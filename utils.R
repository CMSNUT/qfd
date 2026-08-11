#==============================================================================
# Clean CID: keep only numeric CIDs, remove invalid/empty values
#==============================================================================
clean_cids <- function(cid_vec) {
  cid_vec <- as.character(cid_vec)
  cid_vec <- gsub("[[:space:]]", "", cid_vec)  # Remove all whitespace
  cid_vec <- ifelse(grepl("^[0-9]+$", cid_vec), cid_vec, NA)  # Keep digits only
  return(cid_vec[!is.na(cid_vec)])
}

#==============================================================================
# Batch query PubChem API with retry & exponential backoff
#==============================================================================
fetch_batch <- function(properties, batch_cids, retries = 10, delay = 0.5) {
  # Clean input CIDs
  batch_cids <- clean_cids(batch_cids)
  if (length(batch_cids) == 0) return(data.frame())
  
  # Build API URL
  cid_str <- paste(batch_cids, collapse = ",")
  prop_str <- paste(properties, collapse = ",")
  api_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                    cid_str, "/property/", prop_str, "/JSON")
  
  # Retry loop
  for (attempt in seq_len(retries)) {
    result <- tryCatch({
      Sys.sleep(delay)
      res <- httr::GET(api_url)
      
      # Success
      if (httr::status_code(res) == 200) {
        data <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
        df <- data$PropertyTable$Properties
        if (!is.null(df) && nrow(df) > 0) {
          df$CID <- as.character(df$CID)
          return(df)
        } else {
          return(data.frame())
        }
      }
      
      # Server overload: exponential backoff
      if (httr::status_code(res) == 503) {
        warning(sprintf("503 | CID batch: %s | retry %s", cid_str, attempt))
        Sys.sleep(delay * 2^attempt)
        next
      }
      
      # Other HTTP errors
      warning(sprintf("Request failed | code: %s", httr::status_code(res)))
      return(data.frame())
      
    }, error = function(e) {
      warning(sprintf("Error: %s", e$message))
      return(NULL)
    })
    
    if (!is.null(result) && nrow(result) > 0) return(result)
  }
  
  warning(sprintf("All retries failed | CID batch: %s", cid_str))
  return(data.frame())
}


#==============================================================================
# download SDF 
#==============================================================================

download_all_sdf <- function(cid_list, output_dir = "sdf_files") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  results <- data.frame()
  total <- length(cid_list)
  
  for (i in seq_along(cid_list)) {
    cid <- cid_list[i]
    cat(sprintf("[%d/%d] Downloading CID %s ... ", i, total, cid))
    
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/SDF")
    
    # Try 3D
    res <- httr::GET(url, query = list(record_type = "3d"))
    
    if (status_code(res) == 200) {
      out_file <- file.path(output_dir, paste0(cid, ".sdf"))
      writeBin(httr::content(res, "raw"), out_file)
      cat("OK (3D)\n")
      results <- rbind(results, data.frame(CID = cid, status = "success", type = "3d"))
      next
    }
    
    # Try 2D
    res <- httr::GET(url, query = list(record_type = "2d"))
    if (status_code(res) == 200) {
      out_file <- file.path(output_dir, paste0(cid, ".sdf"))
      writeBin(httr::content(res, "raw"), out_file)
      cat("OK (2D)\n")
      results <- rbind(results, data.frame(CID = cid, status = "success", type = "2d"))
      next
    }
    
    cat("FAILED\n")
    results <- rbind(results, data.frame(CID = cid, status = "failed", type = "none"))
    
    Sys.sleep(0.2)
  }
  
  # Summary
  cat("\n========== Summary ==========\n")
  cat("Total:", nrow(results), "\n")
  cat("Success (3D):", sum(results$type == "3d"), "\n")
  cat("Success (2D):", sum(results$type == "2d"), "\n")
  cat("Failed:", sum(results$status == "failed"), "\n")
  
  write.csv(results, file.path(output_dir, "download_report.csv"), row.names = FALSE)
  
  return(results)
}


#==============================================================================
# Read compound-target files (CSV/XLSX) from folder, add CID from filename
#==============================================================================
read_compound_target <- function(folder, skip_row = 0) {
  if (!dir.exists(folder)) stop("Folder not found: ", folder)
  
  # 1. 列出所有支持的文件
  files <- c(
    list.files(folder, "\\.csv$",  full.names = TRUE, ignore.case = TRUE),
    list.files(folder, "\\.tsv$",  full.names = TRUE, ignore.case = TRUE),
    list.files(folder, "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  )
  
  if (length(files) == 0) {
    warning("No CSV/TSV/XLSX files found.")
    return(data.frame())
  }
  
  # 2. 读取单个文件（自动选择最优解析方式）
  read_one_file <- function(file_path, skip_row) {
    cid <- tools::file_path_sans_ext(basename(file_path))
    ext <- tolower(tools::file_ext(file_path))
    
    df <- NULL
    
    if (ext == "csv") {
      # 优先使用 readr，没有则回退到基础函数
      if (requireNamespace("readr", quietly = TRUE)) {
        df <- readr::read_csv(file_path, skip = skip_row, show_col_types = FALSE)
      } else {
        df <- read.csv(file_path, stringsAsFactors = FALSE, 
                       check.names = FALSE, skip = skip_row, row.names = NULL)
      }
      
    } else if (ext == "tsv") {
      # 优先 data.table (最快) -> 其次 readr -> 最后基础函数
      if (requireNamespace("data.table", quietly = TRUE)) {
        df <- data.table::fread(file_path, skip = skip_row, data.table = FALSE)
      } else if (requireNamespace("readr", quietly = TRUE)) {
        df <- readr::read_tsv(file_path, skip = skip_row, show_col_types = FALSE)
      } else {
        # 基础函数读取 TSV 必须指定 sep 和 header
        df <- read.table(file_path, sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE, check.names = FALSE, 
                         skip = skip_row)
      }
      
    } else if (ext == "xlsx") {
      # Excel 必须依赖 readxl，没装就直接报错
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Package 'readxl' is required for XLSX files. Please install it.")
      }
      df <- as.data.frame(
        readxl::read_excel(file_path, .name_repair = "minimal", skip = skip_row),
        stringsAsFactors = FALSE
      )
      
    } else {
      stop("Unsupported file type: ", ext)
    }
    
    # 统一转为 data.frame（防止 tibble 或 data.table 干扰后续 rbind）
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    # 添加来源列
    df <- cbind(CID = cid, df)
    return(df)
  }
  
  # 3. 批量读取并合并
  all_dfs <- lapply(files, function(f) read_one_file(f, skip_row))
  
  # 安全合并（自动补齐列）
  if (requireNamespace("dplyr", quietly = TRUE)) {
    combined <- dplyr::bind_rows(all_dfs)
  } else {
    all_cols <- unique(unlist(lapply(all_dfs, names)))
    combined <- do.call(rbind, lapply(all_dfs, function(df) {
      missing <- setdiff(all_cols, names(df))
      for (col in missing) df[[col]] <- NA
      df[, all_cols, drop = FALSE]
    }))
  }
  
  return(combined)
}


get_active_aids <- function(cid, retry = 3) {
  for (attempt in 1:retry) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                  cid, "/aids/JSON/?aids_type=active")
    response <- httr::GET(url)
    
    if (status_code(response) == 200) {
      data <- jsonlite::fromJSON(content(response, "text", encoding = "UTF-8"))
      aids <- data$InformationList$Information$AID
      if (is.null(aids)) aids <- NA
      return(aids)
    } else if (status_code(response) == 404) {
      return(NA)
    } else {
      Sys.sleep(1)
    }
  }
  warning(paste("Failed to retrieve data for CID:", cid))
  return(NA)
}
