# Synergistic Mechanisms of *Aconitum carmichaelii* and *Astragalus membranaceus* Combination Against Chronic Heart Failure Due to Dilated Cardiomyopathy: Insights from Multi-Omics and Causal Inference

## Overview
This open repository contains fully reproducible analytical R scripts, standardized processed datasets, detailed parameter records, publication-ready visualization pipelines, and complete methodological documentation for the study titled *Synergistic Mechanisms of *Aconitum carmichaelii* and *Astragalus membranaceus* Combination Against Chronic Heart Failure Due to Dilated Cardiomyopathy: Insights from Multi-Omics and Causal Inference*.

The core innovation of this research integrates multi-layer omics profiling with genetic causal inference, establishing a translational bridge between modern bioinformatic evidence and traditional Chinese medicine (TCM) meridian tropism theory. We systematically identify 24 blood-absorbed bioactive substances, characterize organ-specific target modules, dissect cell-type-specific regulatory signatures, and unravel synergistic anti-heart-failure signaling pathways of the medicinal plant combination consisting of *Aconitum carmichaelii* (Fuzi) and *Astragalus membranaceus* (Huangqi). All scripts are comprehensively annotated, and every data processing workflow is fully traceable to ensure complete reproducibility of analytical results.

## Core Analytical Workflow Implemented in This Repo
All pipelines are self-contained and executable under R (v4.5.3), covering the following modules:
  1. **Identification  of blood-absorbed active ingredients & drug target prediction**
   Twenty-four in-vivo absorbed components of the *Aconitum carmichaelii*–*Astragalus membranaceus* combination were screened from published HPLC-TOF-MS metabolomics data, classified into diterpene alkaloids derived from *Aconitum carmichaelii* and isoflavones derived from *Astragalus membranaceus*. Potential protein targets of all ingredients were predicted via three authoritative databases: SwissTargetPrediction, PharmMapper and SEA.
2. **Disease target collection & intersection screening**
   Integrated DCM-CHF disease targets were compiled from multiple sources: differentially expressed genes from bulk transcriptome dataset GSE57338, WGCNA hub genes, risk loci from HERMES2 DCM GWAS, and experimentally validated pathogenic genes retrieved from peer-reviewed literatures. Intersection targets shared by herbal ingredients and DCM were extracted for subsequent downstream analyses.
3. **Multi-organ target localization & TCM meridian tropism validation**
   Tissue-wide expression profiling of core targets was performed using BioGPS and Human Protein Atlas (HPA). All causal core targets were stratified into two independent functional modules corresponding to the core TCM efficacies of *Astragalus membranaceus* (Qi-tonifying) and *Aconitum carmichaelii* (Yang-warming), which perfectly match their classic TCM meridian tropism characteristics.
4. **GO & KEGG functional enrichment analysis**
   Gene set enrichment analyses covering Gene Ontology (biological process, cellular component, molecular function) and KEGG signaling pathways were conducted based on drug-disease intersecting targets.
5. **Immune infiltration quantification**
   Computational deconvolution algorithms were applied to calculate relative abundance of diverse infiltrated immune cell subpopulations in myocardial tissues of DCM patients.
6. **Genetic causal inference pipeline**
   Two-sample Mendelian randomization (MR) combined with Bayesian colocalization analysis was implemented. GTEx v8 left ventricular eQTL data and HERMES2 DCM GWAS summary statistics were adopted to screen causal hub genes and exclude confounding irrelevant correlative signals.
7. **Single-cell RNA-seq analysis & ssGSEA scoring**
   Public single-cell dataset GSE183852 was utilized to characterize cell-type-specific target expression patterns; single-sample gene set enrichment analysis (ssGSEA) quantified overall pathway activity across individual samples.
8. **Molecular docking validation**
   Online tool CB-DOCK2 was used for molecular docking simulation, with results invoked and processed via R system commands to quantify binding affinity between representative bioactive ingredients and causal core targets.
9. **Batch generation of publication figures and supplementary tables**
   Standardized plotting and table export scripts were developed, fully compliant with MDPI journal formatting standards for manuscript visualization.

## Repository Directory Structure
```
qfd/
├── data/                     
│   ├── raw                   # Links of public datasets, index of raw GEO/GWAS/eQTL files
│   ├── processed             # Cleaned gene expression matrices, target lists, SMILES of ingredients
├── scripts/
│   ├── 00_required_packages.qmd  # One-click installation script for all dependent R packages
│   ├── 01_ingredient_identification_target_prediction.R
│   ├── 02_bulk_rna_deg_wgcna.R
│   ├── 03_common_targets_hub_genes.R
│   ├── 04_organ_localization_meridian_tropism.R
│   ├── 05_functional_enrichment.R
│   ├── 06_immune_infiltration_quantification.R
│   ├── 07_mr_colocalization.R
│   ├── 08_single_cell_analysis_ssGSEA.R
│   ├── 09_molecular_docking_analysis.R
│   └── 10_figures_and_tables.R
├── figures/                  # Exported main publication figures (PDF / PNG format)
├── supplementary/
│   ├── figures/              # Supplementary figures for manuscript
│   └── tables/               # Supplementary result tables (CSV / XLSX)
└── README.md
```

## Quick Start Guide
1. Clone this repository to your local machine
```bash
git clone https://github.com/cmsnut/qfd.git
cd qfd
```

## R Package Installation Code
- Recommended R version: 4.5.3. All analytical pipelines in this repository are developed purely in R, no virtual environment required.
- Execute the unified installation script to batch install all required packages: `scripts/00_required_packages.qmd`

## Data and External Tools
### Public Online Databases
- GEO (https://www.ncbi.nlm.nih.gov/geo/)
- GTEx v8 (https://gtexportal.org/home/)
- HERMES dilated cardiomyopathy 2024 GWAS (European ancestry): Heart Failure Molecular Epidemiology for Therapeutic Targets (HERMES 2.0) GWAS (https://www.kp4cd.org/node/1482)
  GWAS summary statistics for heart failure were sourced from the HERMES 2.0 Consortium, covering 6,001 strict DCM-CHF cases and 449,382 control participants.
- BioGPS (http://biogps.org/)
- HPA (https://www.proteinatlas.org/)
- SwissTargetPrediction (https://www.swisstargetprediction.ch/)
- PharmMapper (https://www.lilab-ecust.cn/pharmmapper/)
- SEA (https://seadev.docking.org/)

### External Docking Tool
- CB-DOCK2 (http://183.56.231.194:8001/cb-dock2/index.php): Online platform for molecular docking simulation

## Standard Analytical Pipeline
1. Download all raw public omics and GWAS datasets following access links recorded in `data/raw`
2. Execute analysis scripts sequentially in numerical order under the `scripts/` folder
3. Main manuscript figures are automatically exported to the `figures/` folder; all supplementary figures and result tables are saved to subfolders within `supplementary/`

## Key Scientific Highlights
1. Integrate bulk transcriptomics and single-cell RNA-seq to systematically dissect tissue-wide and cell-type-specific regulatory mechanisms exerted by the combination of *Aconitum carmichaelii* and *Astragalus membranaceus* against dilated cardiomyopathy-induced chronic heart failure.
2. Combine multi-organ target localization profiling with TCM meridian tropism theory, confirming that *Astragalus membranaceus* and *Aconitum carmichaelii* regulate two distinct functional target modules, consistent with their traditional meridian distribution and clinical therapeutic efficacies.
3. Apply rigorous two-sample Mendelian randomization coupled with Bayesian colocalization causal inference to distinguish causal disease-driving hub genes from irrelevant correlative genes, greatly improving the biological reliability of core pharmacodynamic targets compared with conventional network pharmacology frameworks.
4. Systematically identify and characterize 24 blood-absorbed active ingredients from the *Aconitum carmichaelii*–*Astragalus membranaceus* combination, and verify direct molecular binding interactions between representative bioactive components and causal core targets via CB-DOCK2 molecular docking simulation.
5. Deliver a standardized, fully reproducible multi-omics and causal inference analytical workflow for research on medicinal plant pairs, which can be flexibly adapted for mechanistic investigations of other cardiovascular herbal combinations.

## Citation
If you reuse the analytical scripts, pipeline framework or research logic from this repository, please cite our manuscript upon formal publication:
> Synergistic Mechanisms of *Aconitum carmichaelii* and *Astragalus membranaceus* Combination Against Chronic Heart Failure Due to Dilated Cardiomyopathy: Insights from Multi-Omics and Causal Inference. [Journal Name], Year.

## License
This project is released under the **MIT License**. Free academic reuse, modification and secondary development are permitted with appropriate citation of this repository and associated manuscript. Commercial application requires prior communication with the corresponding author.

## Contact
For script bugs, dataset access issues or academic collaboration proposals, please submit an Issue in this repository or contact the corresponding author via email: [Your Email Address]

## Full Online Documentation
👉 Access complete detailed documentation for all analytical pipelines here: https://CMSNUT.github.io/qfd


