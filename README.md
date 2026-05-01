# Chemical Genomics Analysis

This repository contains the code and tools for chemical genomics fitness analysis of bacterial isolates, including figure generation scripts and a machine-learning fitness predictor.

---

## Contents

| Component | Description |
|---|---|
| `Figures_FitnessValue.R` | R script to generate all fitness, GWAS, AMR, pan-genome, and ML figures |
| `fitness-predictor` | CLI tool to predict bacterial fitness from genome FASTA (conda / Docker) |

---

## Figures

### Data availability

Input data are not stored in this repository due to file size. Download from Zenodo:

**Data DOI:** https://doi.org/10.5281/zenodo.19931379

After downloading, place the `Data/` folder in the project root:

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

Unzip the compressed folders before running the script:

```bash
cd Data
unzip assay_results_full.zip
unzip results_prediction.zip
unzip panGenome.zip
unzip SNP_GWAS.zip
cd ..
```

### Data files

| Path | Description |
|---|---|
| `absolute_fitness_ML.csv` | Main absolute-fitness matrix for ML and condition-level analyses |
| `AMR_genes.csv` | AMR gene matrix used for resistance feature analysis and plotting |
| `AMRFinder_results.tsv` | AMRFinder output with detected AMR, virulence, and related gene annotations |
| `condition_group.csv` | Metadata linking conditions to categories, chemicals, and chemical classes |
| `gene_1.csv` | Isolate/gene annotation table for pan-genome and prediction analyses |
| `gene_presence_absence.Rtab` | Roary-style gene presence/absence matrix |
| `group_3530_presence_absence.tsv` | Group-level gene presence/absence table |
| `GWAS_results_SNPs_Pangenome.csv` | Combined SNP and pan-genome GWAS hit table |
| `ML_clusters.csv` | ML cluster assignments for cluster-level analysis |
| `pangenome_pliotropy.csv` | Genes associated with multiple conditions or traits |
| `pan_genome_reference.fa` | Pan-genome FASTA reference |
| `predictive_totall_importance_1.csv` | Aggregated feature-importance table |
| `SHAP_important_features.csv` | SHAP-based feature importance from ML models |
| `SNP_GWAS_annot.csv` | Annotation table for SNP GWAS hits |
| `summary_ARG.tsv` | ARG/plasmid annotation summary |
| `total_distance_matrix.csv` | Pairwise distance matrix for clustering and visualisation |
| `tot_results.tsv` | Combined results table for downstream analysis |
| `conditions/` | Plate-layout files mapping colony positions to isolate identifiers |
| `assay_results_full.zip` | Raw IRIS colony assay result files |
| `results_prediction.zip` | ML prediction result files |
| `panGenome.zip` | Pan-genome GWAS result files |
| `SNP_GWAS.zip` | SNP GWAS result files |

### Required R packages

```r
install.packages(c(
  "tidyverse", "ggplot2", "dplyr", "tidyr", "tibble", "stringr",
  "scales", "ggpointdensity", "diptest", "Rtsne", "mclust",
  "corrplot", "glmnet", "pheatmap", "lme4", "lmerTest",
  "vegan", "ape", "phangorn", "randomcoloR"
))

install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "KEGGREST", "ggtree", "treeio"))
```

### How to run

```bash
# 1. Clone the repository
git clone https://github.com/Sara-Iftikhar/ChemGenomics.git
cd ChemGenomics

# 2. Download data from Zenodo and place Data/ in the project root

# 3. Unzip compressed data
cd Data && unzip assay_results_full.zip && unzip results_prediction.zip
unzip panGenome.zip && unzip SNP_GWAS.zip && cd ..

# 4. Run the script
Rscript Figures_FitnessValue.R
```

Or open `Figures_FitnessValue.R` in RStudio and run section by section.

### Path note

The script uses relative paths (`./Data/`). If absolute paths like `/Data/` are used, mount accordingly:

```bash
docker run -v "$PWD/Data:/Data" <image_name>
```

---

## Fitness Predictor

Predict the absolute fitness of a bacterial strain under different growth conditions from a genome FASTA file using NGBoost probabilistic regression.

### File structure

```
fitness-predictor/
├── predict.py                  # Main prediction script
├── environment.yml             # Conda environment definition
├── requirements.txt            # Pip dependencies
├── pan_genome_reference.fa     # BLAST reference (pan-genome genes)
├── gene_annotations.csv        # Gene annotations (used by --explain)
└── models/                     # 160 trained NGBoost models (.joblib)
```

### Installation

**Option 1 — Conda**

```bash
conda install -c gzhoubioinf -c bioconda -c conda-forge fitness-predictor
```

Then download models and reference data (one time):

```bash
fitness-predict --download-data ~/fitness_data
cd ~/fitness_data
```

Run predictions from `~/fitness_data`:

```bash
# List all available conditions (160 total)
fitness-predict --list

# Predict fitness for a strain
fitness-predict -c <condition> -p /path/to/genome.fasta

# Predict with gzipped FASTA
fitness-predict -c <condition> -p /path/to/genome.fasta.gz

# Custom BLAST thresholds (default: identity=95, coverage=80)
fitness-predict -c <condition> -p /path/to/genome.fasta --identity 80 --coverage 60

# Explain prediction with SHAP gene contributions
fitness-predict -c <condition> -p /path/to/genome.fasta --explain

# Show top 30 genes instead of the default 20
fitness-predict -c <condition> -p /path/to/genome.fasta --explain --top-n 30
```

**Option 2 — Docker**

```bash
docker pull gzhoubioinf09/fitness-predictor:v0.2.0

# List available conditions
docker run --rm gzhoubioinf09/fitness-predictor:v0.2.0 --list

# Predict fitness (e.g. condition=Doripenem0125ugml, FASTA in current directory)
docker run --rm -v $(pwd):/data gzhoubioinf09/fitness-predictor:v0.2.0 -c <condition> -p /data/genome.fasta

# Custom BLAST thresholds
docker run --rm -v $(pwd):/data gzhoubioinf09/fitness-predictor:v0.2.0 -c <condition> -p /data/genome.fasta --identity 80 --coverage 60

# SHAP explanation
docker run --rm -v $(pwd):/data gzhoubioinf09/fitness-predictor:v0.2.0 -c <condition> -p /data/genome.fasta --explain

# SHAP explanation with top 30 genes
docker run --rm -v $(pwd):/data gzhoubioinf09/fitness-predictor:v0.2.0 -c <condition> -p /data/genome.fasta --explain --top-n 30
```

### Output

```
Condition      : Doripenem0125ugml
FASTA          : genome.fasta
Genes detected : 412 / 4751
Fitness (mean) : 1.0842
95% CI         : [0.8901, 1.2783]

Result saved   : genome_Doripenem0125ugml_prediction.csv
```

- **Fitness > 1** — strain grows better than the median strain on the plate
- **Fitness ≈ 1** — strain grows at the median rate
- **Fitness < 1** — strain grows worse than the median strain on the plate

Two CSV files are saved:

| File | Contents |
|---|---|
| `genome_<condition>_prediction.csv` | One-row summary: condition, fasta, genes_hit, total_genes, fitness_mean, ci95_lower, ci95_upper |
| `genome_<condition>_shap.csv` | One row per gene: gene, shap, effect, annotation *(only with `--explain`)* |

With `--explain`, SHAP gene contributions are also printed:

```
SHAP (baseline fitness with no genes: 0.7485)
  Gene                           SHAP      Effect      Annotation
  ---------------------------- ----------  ----------  ----------
  blaCTXM~~~blaCTXM15          +0.14676  promotes    CTX-M family extended-spectrum class A beta-lactamase
  group_5082                   -0.08161  reduces     Transposase
  group_8522                   +0.07659  promotes    ...
```

- **Baseline** — predicted fitness of a hypothetical strain with no pan-genome genes
- **SHAP value** — additive contribution of each present gene; values sum exactly to `fitness − baseline`
- **promotes / reduces** — whether carrying that gene raises or lowers fitness in this condition

### How it works

1. The input genome is BLASTed against a pan-genome reference (`pan_genome_reference.fa`)
2. Gene presence/absence is determined from BLAST hits (identity ≥ 95%, coverage ≥ 80% by default)
3. A pre-trained NGBoost model predicts fitness and a 95% confidence interval from the gene vector

Models were trained on 160 conditions using GWAS-significant genes (Benjamini-Hochberg p < 0.05) as features and Spearman rho as the validation metric (80/20 train/test split).

### Accepted input formats

| Extension | Supported |
|---|---|
| `.fasta`, `.fa`, `.fna` | Yes |
| `.fasta.gz`, `.fa.gz`, `.fna.gz` | Yes (decompressed automatically) |
| `.fasta.bz2`, `.fasta.zst` | No |

---

## Contact

This is a joint project developed by the **Infectious Disease Epidemiology Lab (KAUST)** and the **Banzhaf Lab (Newcastle University)**.

| Name | Role | Email |
|---|---|---|
| Danesh Moradigaravand | PI, KAUST | danesh.moradigaravand@kaust.edu.sa |
| Manuel Banzhaf | PI, Newcastle University | manuel.banzhaf@newcastle.ac.uk |
