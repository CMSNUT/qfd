if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("progressr")) install.packages("progressr")
if (!require("PubChemR")) install.packages("PubChemR")  

library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(progressr)
library(PubChemR)

#' get_active_aids 函数
#'
#' 获取单个 CID 的活性 AID 列表
#'
#' @param cid 整数或字符，PubChem CID
#' @param retry 重试次数
#' @param delay 重试间隔（秒）
#' @return 整数向量
get_active_aids <- function(cid, retry = 1, delay = 1) {
  # 参数校验
  stopifnot(length(cid) == 1)
  cid_char <- as.character(cid)
  
  for (attempt in seq_len(retry)) {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
      cid_char,
      "/aids/JSON?aids_type=active"
    )
    resp <- tryCatch(GET(url), error = function(e) NULL)
    
    if (!is.null(resp) && status_code(resp) == 200) {
      data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
      aids <- data$InformationList$Information$AID
      if (is.null(aids)) return(integer(0))
      # 确保是向量（有时是列表）
      if (is.list(aids)) aids <- unlist(aids)
      return(as.integer(aids))
    } else if (!is.null(resp) && status_code(resp) == 404) {
      return(integer(0))
    }
    Sys.sleep(delay)
  }
  # warning("Failed to retrieve active AIDs for CID: ", cid_char)
  return(integer(0))
}


#' 批量获取多个 CID 的活性 AID，返回数据框
#'
#' @param cids 整数或字符向量，PubChem CID 列表
#' @param retry 每个 CID 请求的重试次数
#' @param delay 每次请求后的等待时间（秒）
#' @param verbose 是否显示进度信息（默认 TRUE）
#' @return 数据框，包含列 CID 和 AID（长格式，每个 CID-AID 一行）
#' @examples
#' \dontrun{
#' cids <- c(2244, 3672, 528882)
#' df <- get_cid_aid_df(cids)
#' print(df)
#' }
get_cid_aid_df <- function(cids, retry = 1, delay = 1, verbose = TRUE) {
  # 确保加载 progressr
  if (!requireNamespace("progressr", quietly = TRUE)) {
    stop("Package 'progressr' is required.")
  }
  
  # 去重并转为字符
  cids_unique <- unique(as.character(cids))
  total <- length(cids_unique)
  if (total == 0) {
    warning("No valid CIDs provided.")
    return(data.frame(CID = character(), AID = integer(), stringsAsFactors = FALSE))
  }
  
  if (verbose) message("Processing ", total, " CIDs...")
  
  # 存储结果的列表
  result_list <- list()
  
  # 使用 progressr 显示进度
  with_progress({
    p <- progressor(along = cids_unique)
    
    for (i in seq_along(cids_unique)) {
      cid <- cids_unique[i]
      aids <- get_active_aids(cid, retry = retry, delay = delay)
      
      # 如果有 AID，构建临时数据框
      if (length(aids) > 0) {
        temp_df <- data.frame(
          CID = cid,
          AID = aids,
          stringsAsFactors = FALSE
        )
        result_list[[i]] <- temp_df
      } else {
        result_list[[i]] <- NULL
      }
      
      p(sprintf("CID=%s", cid))
      Sys.sleep(delay)  
    }
  })
  
  # 合并所有结果
  final_df <- bind_rows(result_list)

  if (verbose) message("Done. Retrieved ", nrow(final_df), " CID-AID pairs.")
  return(final_df)
}


get_assays_df <- function(cid, retry = 1, delay = 1) {
  # 参数校验
  stopifnot(length(cid) == 1)
  cid_char <- as.character(cid)
  
  for (attempt in seq_len(retry)) {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
      cid_char,
      "/assaysummary/JSON"
    )
    resp <- tryCatch(GET(url), error = function(e) NULL)
    
    if (!is.null(resp) && status_code(resp) == 200) {
      data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
      if (length(data$Table$Row$Cell) >0) {
        df <- as.data.frame(do.call(rbind, data$Table$Row$Cell))
        colnames(df) <- data$Table$Columns$Column
      } else {
        df <- NULL
      }
    } else if (!is.null(resp) && status_code(resp) == 404) {
      df <- NULL
    }
    Sys.sleep(delay)
  }
  
  return(df)
}


# 
get_assays_by_cids <- function(cids, retry = 1, delay = 2, verbose = TRUE) {
  # 确保加载 progressr
  if (!requireNamespace("progressr", quietly = TRUE)) {
    stop("Package 'progressr' is required.")
  }
  
  # 去重并转为字符
  cids_unique <- unique(as.character(cids))
  total <- length(cids_unique)
  if (total == 0) {
    warning("No valid CIDs provided.")
    return(NULL)
  }
  
  if (verbose) message("Processing ", total, " CIDs...")
  
  # 存储结果的列表
  result_list <- list()
  
  # 使用 progressr 显示进度
  with_progress({
    p <- progressor(along = cids_unique)
    
    for (i in seq_along(cids_unique)) {
      cid <- cids_unique[i]
      df <- get_assays_df(cid, retry = retry, delay = delay)
      if (!is.null(df)) {
        df <- df|>
          mutate(across(everything(), as.character))
        result_list[[i]] <- df
      } else {
        result_list[[i]] <- NULL
      }

      p(sprintf("CID=%s", cid))
      Sys.sleep(delay)
    }
  })
  
  # 合并所有结果
  final_df <- bind_rows(result_list)
  
  if (verbose) message("Done. Retrieved ", nrow(final_df), " CID-AID pairs.")
  return(final_df)
}


get_cid_targets <- function(cids, retry = 1, delay = 1, verbose = TRUE) {
  # 确保加载 progressr
  if (!requireNamespace("progressr", quietly = TRUE)) {
    stop("Package 'progressr' is required.")
  }
  
  # 去重并转为字符
  cids_unique <- unique(as.character(cids))
  total <- length(cids_unique)
  if (total == 0) {
    warning("No valid CIDs provided.")
    # return(data.frame(CID = character(), AID = integer(), stringsAsFactors = FALSE))
  }
  
  handlers("txtprogressbar")  # 或 "cli" 更美观
  
  if (verbose) message("Processing ", total, " CIDs...")
  
  # 存储结果的列表
  result_list <- list()
  
  # 使用 progressr 显示进度
  with_progress({
    p <- progressor(along = cids_unique)
    
    for (i in seq_along(cids_unique)) {
      cid <- cids_unique[i]
      
      query_list <- list(
        download = "*",
        collection = "consolidatedcompoundtarget",
        order = list("cid,asc"),
        start = 1,
        limit = 10000000,
        downloadfilename = paste0("pubchem_cid_",cid,"_consolidatedcompoundtarget"),
        where = list(
          ands = list(
            list(cid=cid)
          )
        )
      )
      
      query_json <- toJSON(query_list, auto_unbox = TRUE, pretty = FALSE)
      query_encoded <- URLencode(query_json, reserved = TRUE)
      base_url <- "https://pubchem.ncbi.nlm.nih.gov/sdq/sphinxql.cgi"
      full_url <- paste0(base_url, "?infmt=json&outfmt=json&query=", query_encoded)
      
      
      # 重试循环 
      success <- FALSE
      last_error <- NULL
      
      for (attempt in 0:retry) {
        tryCatch({
          # 发送请求
          response <- GET(full_url, add_headers(`User-Agent` = "Mozilla/5.0"))
          
          # 检查状态码
          if (status_code(response) != 200) {
            stop(paste("HTTP状态码:", status_code(response)))
          }
          
          # 解析 JSON
          text_raw <- content(response, "text", encoding = "UTF-8")
          
          # 尝试直接解析
          data <- tryCatch({
            fromJSON(text_raw, flatten = TRUE)
          }, error = function(e) {
            # 解析失败，修复格式后再试
            fixed_text <- gsub("\\}\\s*\\{", "},{", text_raw)
            fromJSON(fixed_text, flatten = TRUE)
          })
          
          # 成功
          if (length(data) >0) {
            message(cid, " 数据获取成功！共 ", nrow(data), " 条记录。")
            result_list[[i]] <- data |>
              mutate(across(everything(), as.character))
          } else {
            message(cid, " 数据获取成功！共 0 条记录。")
            result_list[[i]] <- NULL
          }
          
          success <- TRUE
          break   # 成功则跳出重试循环
          
        }, error = function(e) {
          # 记录错误信息
          last_error <<- paste("尝试", attempt + 1, "失败:", e$message)
          message(cid, " ", last_error)
          
          # 如果还有重试次数，等待后继续
          if (attempt < max_retries) {
            wait_time <- retry_delay * (2 ^ attempt)  # 指数退避：2, 4, 8...
            message("等待 ", wait_time, " 秒后重试...")
            Sys.sleep(wait_time)
          } else {
            # 已达最大重试次数
            message(cid, " 已达到最大重试次数，放弃。")
          }
        })
        
        # 如果成功，跳出重试循环（已被 break 处理）
      }
      
      # 如果最终仍未成功，记录错误并置空
      if (!success) {
        result_list[[i]] <- NULL
        error_log[[i]] <- last_error
      }
      
      
      p(sprintf("CID=%s", cid))
      Sys.sleep(delay)  
    }
  })
  
 
  
  # 合并所有结果
  final_df <- bind_rows(result_list)

  if (verbose) message("Done. Retrieved ", nrow(final_df), " CID-Target pairs.")
  return(final_df)
}


get_gene_symbol_by_aid <- function(aid, delay = 1) {
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/",
    as.character(aid),
    "/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/JSON")
  
  resp <- tryCatch(GET(url), error = function(e) NULL)
  
  if (!is.null(resp) && status_code(resp) == 200) {
    data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
    aid <- data$InformationList$Information$AID
    if (is.null(aids)) return(integer(0))
    # 确保是向量（有时是列表）
    if (is.list(aids)) aids <- unlist(aids)
    return(as.integer(aids))
  } else if (!is.null(resp) && status_code(resp) == 404) {
    return(integer(0))
  }
  Sys.sleep(delay)
}
  

get_gene_symbol_by_geneids <- function(geneids, delay = 1, verbose = TRUE) {
  if (!requireNamespace("progressr", quietly = TRUE)) {
    stop("Package 'progressr' is required.")
  }
  
  geneids_unique <- unique(as.character(geneids))
  geneids_unique <- Filter(nzchar, geneids_unique)
  total <- length(geneids_unique)
  
  if (total == 0) {
    warning("No valid CIDs provided.")
    return(data.frame(GeneID = character(), 
                      Symbol = character(), 
                      TaxonomyID = character(),
                      stringsAsFactors = FALSE)
          )
  }
  
  if (verbose) message("Processing ", total, " GeneIDs...")
  
  result_list <- list()
  
  with_progress({
    p <- progressor(along = geneids_unique)
    
    for (i in seq_along(geneids_unique)) {
      geneid <- geneids_unique[i]
      url <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/",
        as.character(geneid),
        "/summary/JSON")
      resp <- tryCatch(GET(url), error = function(e) NULL)
  
      if (!is.null(resp) && status_code(resp) == 200) {
        data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
        d <- data$GeneSummaries$GeneSummary
        if (is.null(d)) {
          result_list[[i]] <- NULL
        } else {
          d <- d |> 
            mutate(across(everything(), as.character),
                   across(everything(), trimws)) |>
            dplyr::select(GeneID,Symbol,TaxonomyID)
          result_list[[i]] <- d
        }
        
      } else if (!is.null(resp) && status_code(resp) == 404) {
        result_list[[i]] <- NULL
      }
    p(sprintf("GeneID=%s", geneid))
    Sys.sleep(delay)
    }
  })
  
  # 合并所有结果
  final_df <- bind_rows(result_list)
  
  if (verbose) message("Done. Retrieved ", nrow(final_df), " CID-AID pairs.")
  return(final_df)
}


# library(rentrez)

# GI → Gene Symbol
accession_to_symbol <- function(accession) {
  if (substr(accession, 3, 3) == "_") {
    query = paste0("xref:",accession)
  } else if (grepl("^[0-9]$", substr(accession, 1, 1))) {
    query = paste0("pdb-",substr(accession, 1, 4))
  } else {
    query = accession
  }
  
  resp <- GET(url = "https://rest.uniprot.org/uniprotkb/search",
              query = list( query = query, format = "json"))
  
  # data <- content(resp, "parsed", type = "application/json")
  
  data <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
  
  if (length(data$results) == 0) {
    return(NA)
  }
  
  if (!("genes" %in% names(data$results))) return(NA)
  
  genes_df <- data$results[, c("primaryAccession", "genes")]
  genes_df <- unnest(genes_df, cols = c(genes))
  
  d <- data.frame(
    Accession = rep(accession, nrow(genes_df)),
    Symbol = genes_df$geneName.value,
    TaxonomyID = as.character(data$results$organism.taxonId)
  )
  
  return(d)
}

get_gene_symbol_by_accessions <- function(accessions, delay = 1, verbose = TRUE) {
  
  if (!requireNamespace("progressr", quietly = TRUE)) {
    stop("Package 'progressr' is required.")
  }
  
  accessions_unique <- unique(as.character(accessions))
  accessions_unique <- Filter(nzchar, accessions_unique)
  total <- length(accessions_unique)
 
  if (total == 0) {
    warning("No valid Accessions provided.")
    return(data.frame(Accession= character(), 
                      Symbol = character(), 
                      TaxonomyID = character(),
                      stringsAsFactors = FALSE)
    )
  }
  
  if (verbose) message("Processing ", total, " Accessions...")
  
  result_list <- list()
  
  with_progress({
    p <- progressor(along = accessions_unique)
    
    for (i in seq_along(accessions_unique)) {
      accession <- accessions_unique[i]
      
      df <- accession_to_symbol(accession)
      
      if (is.null(names(df))) {
        result_list[[i]] <- NULL
      } else {
        result_list[[i]] <- df
      }
      
      p(sprintf("Accession=%s", accession))
      Sys.sleep(delay)
    }
    
  })
  
  # 合并所有结果
  final_df <- bind_rows(result_list)
  
  if (verbose) message("Done. Retrieved ", nrow(final_df), " Accession-Symbol pairs.")
  return(final_df)
}


  