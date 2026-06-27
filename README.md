# Integrating Multi-Omics & Causal Inference: Synergistic Mechanisms of Fuzi-Huangqi Herb Pair against DCM-Induced Chronic Heart Failure

## Repository Overview
This open repository stores fully reproducible analytical R scripts, standardized processed datasets, detailed parameter records, publication-grade visualization pipelines, and complete methodological documentation for the manuscript titled *Integrating Multi-Omics and Causal Inference Reveals the Synergistic Mechanisms of the Fuzi-Huangqi Herb Pair Against Chronic Heart Failure Due to Dilated Cardiomyopathy*.

The core innovation of this research integrates multi-layer omics profiling with genetic causal inference, constructing a translational bridge between modern bioinformatic evidence and traditional Chinese medicine (TCM) meridian tropism theory. We systematically characterize blood-absorbed active substances, organ-specific target modules, cell-type regulatory patterns, and synergistic anti-heart-failure signaling pathways of the Fuzi-Huangqi herb pair. All scripts are well-annotated, and every data processing workflow is fully traceable to guarantee reproducible analytical outcomes.

## Core Analytical Workflow Implemented in This Repo
All pipelines are self-contained and executable under R (v4.5.3) environment, included:
1. **Curation of blood-absorbed active ingredients & drug target prediction**
   Twenty-four in-vivo absorbed components of the Fuzi-Huangqi herb pair are screened from published HPLC-TOF-MS metabolomics data, classified into Fuzi-derived diterpene alkaloids and Huangqi-derived isoflavones. Potential protein targets of all ingredients are predicted via SwissTargetPrediction, PharmMapper and SEA.
2. **Disease target collection & intersection screening**
   Integrated DCM-CHF disease targets consist of bulk transcriptome (GSE57338) differentially expressed genes, WGCNA hub genes, HERMES2 DCM GWAS risk loci, and literature-supported pathogenic genes. Common targets shared by herbal ingredients and DCM are extracted for downstream analysis.
3. **Multi-organ target localization & TCM meridian tropism validation**
   Tissue expression profiling of core targets via BioGPS and Human Protein Atlas (HPA). All causal core targets are classified into two functional modules corresponding to Huangqi (Qi-tonifying) and Fuzi (Yang-warming), which match their canonical TCM meridian tropism features.
4. **GO & KEGG functional enrichment analysis**
   Enrichment of biological processes, cellular components, molecular functions and signaling pathways based on intersecting drug-disease targets.
5. **Immune infiltration quantification**
   Deconvolution algorithm to calculate relative proportions of infiltrated immune cells in DCM myocardial tissue.
6. **Genetic causal inference pipeline**
   Two-sample Mendelian randomization (MR) paired with Bayesian colocalization analysis. GTEx v8 left ventricular eQTL data and HERMES2 DCM GWAS summary statistics are utilized to filter causal hub genes and eliminate confounding correlative signals.
7. **Single-cell RNA-seq analysis & ssGSEA scoring**
   Single-cell dataset GSE183852 is used to explore cell-type-specific target expression; ssGSEA quantifies pathway activity across samples.
8. **Molecular docking validation**
   Standalone AutoDock Vina binary invoked via R system commands to calculate binding affinity between representative bioactive ingredients and causal core targets.
9. **Batch generation of publication figures and supplementary tables**
   Standardized plotting and table export scripts compliant with MDPI journal formatting requirements.

## Repository Directory Structure
```
Fuzi-Huangqi-DCM-CHF-MultiOmics-CausalInference/
├── data/                     
│   ├── raw                   # Links of public datasets, index of raw GEO/GWAS/eQTL files
│   ├── processed             # Cleaned gene expression matrices, target lists, SMILES of ingredients
├── scripts/
│   ├── 01_ingredient_target_prediction.R
│   ├── 02_bulk_rna_deg_wgcna.R
│   ├── 03_common_targets_hub_genes.R
│   ├── 04_organ_localization_meridian_tropism.R
│   ├── 05_functional_enrichment.R
│   ├── 06_immune_enrichment.R
│   ├── 07_mr_bayesian_colocalization.R
│   ├── 08_single_cell_analysis_ssGSEA.R
│   ├── 09_molecular_docking_analysis.R
│   └── 10_figures_and_tables.R
├── figures/                  # Exported main publication figures (PDF / PNG format)
├── supplementary/
│   ├── figures/              # Supplementary figures for manuscript
│   └── tables/               # Supplementary result tables (CSV / XLSX)
└── README.md
```

## Software & Environment Setup
All analytical pipelines in this repository are purely built with R language (v4.2.1), no Python scripts are integrated. **No virtual environment is required**. Users can install all dependent packages directly into the global R library.
- Recommended R version: 4.5.3
- Core required packages: WGCNA, clusterProfiler, org.Hs.eg.db, GSVA, immunedeconv, coloc, TwoSampleMR, ggplot2, dplyr, tidyr
- External online tool **CB-DOCK2** (http://183.56.231.194:8001/cb-dock2/index.php) for molecular docking.
- Public online databases adopted: 
  - GEO (https://www.ncbi.nlm.nih.gov/geo/), 
  - GTEx v8 (https://gtexportal.org/home/), 
  - HERMES dilated cardiomyopathy 2024 GWAS: European ancestry, Heart Failure Molecular Epidemiology for Therapeutic Targets (HERMES 2.0) GWAS (https://www.kp4cd.org/node/1482),
  GWAS summary statistics for HF came from the Heart Failure Molecular Epidemiology for Therapeutic Targets (HERMES 2.0) Consortium, including  6,001 cases (strict definition) and 449,382 controls.
  - BioGPS (http://biogps.org/), 
  - HPA (https://www.proteinatlas.org/), 
  - SwissTargetPrediction (https://www.swisstargetprediction.ch/), 
  - PharmMapper (https://www.lilab-ecust.cn/pharmmapper/), 
  - SEA(https://seadev.docking.org/).

### One-Click R Package Installation Code
Run the code snippet below directly in the R console to batch install all required dependencies:
```r
# Install CRAN packages
install.packages(c("dplyr","tidyr","ggplot2","WGCNA","TwoSampleMR","coloc","immunedeconv"))

# Install Bioconductor omics and annotation packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("org.Hs.eg.db", "clusterProfiler", "GSVA"))
```

## Quick Start Guide
1. Clone this repository to your local machine
```bash
git clone https://github.com/[YourUsername]/Fuzi-Huangqi-DCM-CHF-MultiOmics-CausalInference.git
cd Fuzi-Huangqi-DCM-CHF-MultiOmics-CausalInference
```
2. Execute the one-click R installation script above to finish dependency deployment
3. Download all raw public omics and GWAS datasets following links recorded in `data/raw`
4. Run analysis scripts sequentially according to the numeric prefix under the `scripts/` folder
5. Main figures will be saved to `figures/`, while supplementary charts and tables will be output to `supplementary/`

## Key Scientific Highlights
1. Integrate bulk and single-cell transcriptomics to systematically resolve tissue-wide and cell-type-specific regulatory patterns of the Fuzi-Huangqi herb pair against dilated cardiomyopathy-induced chronic heart failure.
2. Combine multi-organ target localization with TCM meridian tropism theory, verifying that Huangqi and Fuzi modulate distinct target modules consistent with their traditional meridian distribution and clinical therapeutic efficacies.
3. Adopt rigorous Mendelian randomization plus Bayesian colocalization causal inference to distinguish disease-causal driving targets from irrelevant correlative genes, significantly improving the reliability of core pharmacodynamic targets compared with conventional network pharmacology workflows.
4. Systematically characterize 24 blood-absorbed active ingredients of the Fuzi-Huangqi herb pair, and validate direct molecular binding interactions between representative bioactive components and causal hub targets via molecular docking.
5. Provide a standardized, fully reproducible multi-omics plus causal inference analytical pipeline for TCM herb pair research, which can be adapted to mechanistic studies of other cardiovascular herbal combinations.

## Citation
If you reuse the scripts, analytical pipeline or research framework from this repository, please cite our manuscript after official publication:
> Integrating Multi-Omics and Causal Inference Reveals the Synergistic Mechanisms of the Fuzi-Huangqi Herb Pair Against Chronic Heart Failure Due to Dilated Cardiomyopathy. [Journal Name], Year.

## License
This project is distributed under the **MIT License**. Free academic reuse, modification and secondary development are permitted with proper citation. Commercial usage requires prior communication with the corresponding author.

## Contact
For script bugs, dataset questions or academic cooperation inquiries, feel free to open an Issue or contact the corresponding author via email: [Your Email Address]

### Modification Notes
1. Fully replaced directory tree with your new structure; removed the obsolete `environment/` folder and related txt file
2. Rewrote the workflow section to strictly match the serial order of your 10 R scripts
3. Split figure outputs into main `figures/` and supplementary subfolders, updated corresponding operation descriptions in Quick Start
4. Removed all Python, virtual environment and renv related content
5. Simplified the package installation code without external txt dependency for lighter use
6. Adjusted folder naming descriptions uniformly to match `data/raw` & `data/processed` standards

## Documentation
👉 Access the complete documentation here: https://CMSNUT.github.io/qfd



