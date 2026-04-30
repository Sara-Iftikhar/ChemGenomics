# Chemical Genomics Fitness Value Analysis

This repository contains the R script required to generate fitness-value, GWAS, AMR, pangenome, and machine-learning related figures for the chemical genomics analysis.

The input data are not stored in this GitHub repository because the data files are large.  
The complete dataset can be downloaded from Zenodo:

**Data DOI:** https://doi.org/10.5281/zenodo.19931379

## Main script

```text
Figures_FitnessValue.R
```

This script performs downstream analyses and plotting using colony fitness measurements, condition metadata, AMR gene results, GWAS outputs, pangenome data, distance matrices, clustering results, and machine-learning prediction results.

## Data availability

All required input data are available on Zenodo:

```text
https://doi.org/10.5281/zenodo.19931379
```

After downloading the data from Zenodo, place the downloaded `Data/` folder in the project root directory so that the folder structure looks like this:

```text
.
├── Figures_FitnessValue.R
├── README.md
└── Data/
    ├── absolute_fitness_ML.csv
    ├── AMR_genes.csv
    ├── AMRFinder_results.tsv
    ├── condition_group.csv
    ├── gene_1.csv
    ├── gene_presence_absence.Rtab
    ├── group_3530_presence_absence.tsv
    ├── GWAS_results_SNPs_Pangenome.csv
    ├── ML_clusters.csv
    ├── pangenome_pliotropy.csv
    ├── pan_genome_reference.fa
    ├── predictive_totall_importance_1.csv
    ├── SHAP_important_features.csv
    ├── SNP_GWAS_annot.csv
    ├── summary_ARG.tsv
    ├── total_distance_matrix.csv
    ├── tot_results.tsv
    ├── conditions/
    │   ├── plate_tr_1.tsv
    │   ├── plate_tr_2.tsv
    │   ├── plate_tr_3.tsv
    │   └── plate_tr_4.tsv
    ├── assay_results_full.zip
    ├── results_prediction.zip
    ├── panGenome.zip
    └── SNP_GWAS.zip
```

## Data files and folders

Each item in the `Data/` folder is described briefly below.

| Path | Description |
|---|---|
| `Data/absolute_fitness_ML.csv` | Main absolute-fitness matrix used for machine-learning and condition-level fitness analyses. |
| `Data/AMR_genes.csv` | Antimicrobial-resistance gene matrix or summary used for AMR feature analysis and plotting. |
| `Data/AMRFinder_results.tsv` | AMRFinder output table containing detected antimicrobial-resistance, virulence, or related gene annotations. |
| `Data/condition_group.csv` | Metadata table linking experimental condition/file names to condition categories, chemicals, and chemical classes. |
| `Data/gene_1.csv` | Isolate or gene annotation table used to link isolates/features with downstream pangenome and prediction analyses. |
| `Data/gene_presence_absence.Rtab` | Roary-style gene presence/absence matrix used as pangenome input for prediction and gene-based analyses. |
| `Data/group_3530_presence_absence.tsv` | Group-level gene presence/absence table used for pangenome or feature-level analysis. |
| `Data/GWAS_results_SNPs_Pangenome.csv` | Combined GWAS hit table containing SNP and pangenome association results used for downstream summaries and figures. |
| `Data/ML_clusters.csv` | Machine-learning cluster assignment table used for cluster-level analysis and plotting. |
| `Data/pangenome_pliotropy.csv` | Pangenome pleiotropy results summarizing genes associated with multiple conditions or traits. |
| `Data/pan_genome_reference.fa` | FASTA reference file for the pan-genome sequences used for annotation or sequence lookup. |
| `Data/predictive_totall_importance_1.csv` | Aggregated predictive feature-importance table used for summary plots of important genes/features. |
| `Data/SHAP_important_features.csv` | SHAP-based feature-importance results from machine-learning models. |
| `Data/SNP_GWAS_annot.csv` | Annotation table for SNP GWAS hits, used to map significant SNPs to genes or functional descriptions. |
| `Data/summary_ARG.tsv` | Summary table of ARG/plasmid-related annotations used to compare resistance profiles and fitness patterns. |
| `Data/total_distance_matrix.csv` | Pairwise distance matrix used for clustering, similarity analysis, or distance-based visualization. |
| `Data/tot_results.tsv` | Total or combined results table used in downstream analysis and plotting. |
| `Data/conditions/` | Folder containing plate-layout files that map colony positions to strain/isolate identifiers for each plate replicate. |
| `Data/conditions/plate_tr_1.tsv` | Plate-layout file for replicate/plate 1. |
| `Data/conditions/plate_tr_2.tsv` | Plate-layout file for replicate/plate 2. |
| `Data/conditions/plate_tr_3.tsv` | Plate-layout file for replicate/plate 3. |
| `Data/conditions/plate_tr_4.tsv` | Plate-layout file for replicate/plate 4. |
| `Data/assay_results_full.zip` | Compressed archive containing raw IRIS colony assay result files for each condition and replicate; unzip before running fitness-processing sections. |
| `Data/results_prediction.zip` | Compressed archive containing machine-learning prediction result files; unzip before running sections that read prediction outputs. |
| `Data/panGenome.zip` | Compressed archive containing pan-genome GWAS result files; unzip before running pangenome GWAS sections. |
| `Data/SNP_GWAS.zip` | Compressed archive containing SNP GWAS result files; unzip before running SNP GWAS plotting sections. |
| `Data/.DS_Store` | macOS system metadata file; not required for the analysis and can be ignored. |

## Compressed folders

Some required inputs are provided as ZIP files in the Zenodo dataset. Unzip them before running the R script.

From the project root, run:

```bash
cd Data
unzip assay_results_full.zip
unzip results_prediction.zip
unzip panGenome.zip
unzip SNP_GWAS.zip
cd ..
```

After unzipping, the data folder should contain the extracted folders/files needed by the script, including assay results, prediction outputs, pan-genome GWAS results, and SNP GWAS results.

## Path note

The folder shown here uses a relative path:

```text
./Data/
```

If `Figures_FitnessValue.R` has been edited to use absolute paths such as:

```r
/Data/...
```

then the script will look for a folder called `/Data` at the root of the system, not inside the current project folder.

To avoid path errors, use one of these options:

### Option 1: Run in an environment where `/Data` points to this project data folder

For example, in a container, mount the folder as `/Data`:

```bash
docker run -v "$PWD/Data:/Data" <image_name>
```

### Option 2: Use a symbolic link if permitted on the system

```bash
ln -s "$PWD/Data" /Data
```

This may require administrator permission depending on the system.

### Option 3: Use relative paths in the R script

If you are running directly from this folder, paths can be adjusted to use:

```r
Data/...
```

instead of:

```r
/Data/...
```

## Required R packages

Install the following R packages before running the script:

```r
install.packages(c(
  "tidyverse",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "scales",
  "ggpointdensity",
  "diptest",
  "Rtsne",
  "mclust",
  "corrplot",
  "glmnet",
  "pheatmap",
  "lme4",
  "lmerTest",
  "vegan",
  "ape",
  "phangorn",
  "randomcoloR"
))
```

Some packages are usually installed through Bioconductor:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler",
  "KEGGREST",
  "ggtree",
  "treeio"
))
```

## How to run

1. Clone this repository:

```bash
git clone https://github.com/Sara-Iftikhar/ChemGenomics.git
cd ChemGenomics
```

2. Download the required data from Zenodo:

```text
https://doi.org/10.5281/zenodo.19931379
```

3. Place the downloaded `Data/` folder inside the project root.

4. Unzip the compressed folders inside `Data/`:

```bash
cd Data
unzip assay_results_full.zip
unzip results_prediction.zip
unzip panGenome.zip
unzip SNP_GWAS.zip
cd ..
```

5. Run the R script from the project root directory:

```bash
Rscript Figures_FitnessValue.R
```

Alternatively, open `Figures_FitnessValue.R` in RStudio and run the script section by section.

```

## Notes

- The data are available from Zenodo and are not included directly in this GitHub repository.
- Unzip all compressed folders before running the relevant sections of the script.
- Make sure the paths in the R script match the actual location of the `Data` folder.
- The script is designed as an analysis/figure-generation workflow and may be easier to run section by section in RStudio.
