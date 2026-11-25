# Chemical Genomics Analysis

This repository contains code and tools for high-throughput chemical genomics screening of *Klebsiella pneumoniae* clinical isolates across 214 growth conditions, including corrected IRIS colony assay files, figure-generation scripts, and an NGBoost-based command-line predictor that estimates strain fitness from genome assemblies.

---

## Contents

| Component | Description |
|---|---|
| `iris_processed_file/` | Manually corrected IRIS colony assay files (`*_modified.iris`) |
| `OriginalAnalysis/Figures_FitnessValue.R` | R script to generate all fitness, GWAS, AMR, pan-genome, and ML figures |
| `fitness-predictor/` | CLI tool to predict bacterial fitness from genome FASTA (conda / Docker) |
| Revision analyses | Scripts and data for the manuscript revision, available on Zenodo (see [Revision Analyses](#revision-analyses)) |

---

## IRIS Processed Files

The `iris_processed_file/` directory contains corrected IRIS colony assay files (`*_modified.iris`). Raw IRIS measurements record zero-size colonies that may represent poor growth, normal growth, over growth rather than true absences. Each zero-size colony was manually reviewed and reclassified before being written back into the corrected `_modified.iris` files used for downstream fitness analysis.

---

## Original Analysis

### Data availability

Input data are not stored in this repository due to file size. Download from Zenodo:

**Data DOI:** https://doi.org/10.5281/zenodo.19931379

`Figures_FitnessValue.R` uses this original release. The script reads its input from the absolute path `/Data/`, so place the downloaded files there (or mount the folder at `/Data`, e.g. `docker run -v "$PWD/Data:/Data" ...`):

```text
/Data/
├── *.csv, *.tsv, *.Rtab, *.fa    # tables listed below
├── conditions/plate_tr_1-4.tsv   # move the four plate_tr_*.tsv files here
└── *.zip                         # unzip before running
```

```bash
cd /Data && mkdir -p conditions && mv plate_tr_*.tsv conditions/
for z in *.zip; do unzip -q "$z"; done
Rscript /path/to/repo/OriginalAnalysis/Figures_FitnessValue.R
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

## Revision Analyses

Analyses added during manuscript revision. Their scripts and data are in the revised Zenodo release (`Revision_*.zip`, `OriginalData/`), in which isolates are identified by **BioSample accession** and only the **1,462 study isolates** and **214 conditions** are kept.

### Setup

Place the revised dataset in `OriginalData/` and extract all `Revision_*.zip` files in the project root:

```text
.
├── OriginalData/
└── Revision/
    ├── Shared_input/          # inputs used by several analyses
    ├── Clonality_analysis/
    ├── EdgeEffect_analysis/
    ├── Opacity_analysis/
    ├── Opacity_ST23/
    ├── Overfitting_analysis/
    ├── Phylogenetic_tree/
    └── Plasmid_analysis/
```

Each analysis folder has one script named after the folder, an `output/` folder and, where needed, an `input/` folder. Run from the project root:

```bash
Rscript Revision/Plasmid_analysis/Plasmid_analysis.R
```

Required R packages: `tidyverse`, `patchwork`, `scales`, `ragg`, `broom`, `fixest`, `progressr`, `ape`, `phangorn`, `ggtree`, `ggtreeExtra` (Bioconductor).

### Analyses

| Folder | Analysis | Main outputs |
|---|---|---|
| `Phylogenetic_tree` | Phylogeny annotated with ST, resistance and virulence | `Figure_S1_*` panels, rooted trees (`.nwk`) |
| `Clonality_analysis` | Genomic vs phenotypic similarity; GWAS robustness to clone pruning | `phenotypic_similarity/`, `clone_pruning_comparison_0.001/` and `_010/` |
| `EdgeEffect_analysis` | Plate edge-effect correction and its effect on GWAS | corrected fitness, `spatial_sensitivity.csv`, `Figure_S15_GWAS_spatial_concordance` |
| `Opacity_analysis` | ST and capsule effects on colony area vs opacity | `ST_capsule_*.csv`, `Figure_Sy_Virulence_Coefficient_Sensitivity` |
| `Opacity_ST23` | Size–opacity concordance and ST23 deviation | `Figure_ST23_Size_Opacity_Concordance`, concordance tables |
| `Overfitting_analysis` | Train–test gaps and feature counts of the fitness models | `FigureS10_Overfitting_FeatureCounts_SD` |
| `Plasmid_analysis` | Replicon vs cargo fitness effects (LMMs, matching, ONT linkage, copy number) | `Figure_S1_LMM` … `Figure_S10_*`, `Table_S*` and `Table_LMM_*` |

`Plasmid_analysis/output/additional/` holds tables produced by code that is not included.

### Input files

Files in each analysis's own `input/` folder (the bold rows list the shared and `OriginalData/` files it also reads):

| File | Description |
|---|---|
| **`Phylogenetic_tree`** | *also reads `OriginalData/tot_results.tsv`* |
| `aln_snps.tree` | FastTree phylogeny of the 1,462 isolates from the core-SNP alignment |
| **`Clonality_analysis`** | *also reads `OriginalData/absolute_fitness_ML.csv`* |
| `pairwise_snp_distances_fixed.tsv` | Pairwise SNP distances between all isolates |
| `merged_gene_gwas_results.csv` | Pan-genome GWAS (Scoary) results per condition, all isolates |
| `merged_scoary_results_0.001.csv`, `merged_scoary_results_010.csv` | The same GWAS after clone pruning at thresholds 0.001 and 0.10 (10 replicates each) |
| **`EdgeEffect_analysis`** | *also reads `OriginalData/assay_results_full.zip` and `plate_tr_1–4.tsv`* |
| `pyseer/original/`, `pyseer/spatial/` | pyseer GWAS per condition on original and spatially corrected fitness |
| **`Opacity_analysis`** | *also reads `Shared_input/Opacity_all_conditions_traits.csv`, `OriginalData/absolute_fitness_ML.csv`* |
| `supp_table_S2.csv` | Supplementary Table S2: Kleborate genotypes per isolate (ST, K/O locus, virulence and resistance scores) |
| **`Opacity_ST23`** | *reads only shared files: `Shared_input/Opacity_all_conditions_traits.csv`, `all_models_mean_sd_five_runs.csv`; `OriginalData/absolute_fitness_ML.csv`, `tot_results.tsv`, `condition_group.csv`* |
| **`Overfitting_analysis`** | *also reads `Shared_input/all_models_mean_sd_five_runs.csv`, `OriginalData/condition_group.csv`* |
| `feature_counts_aggregated.csv` | Number of pan-genome, SNP and unitig features selected per condition (mean ± SD over five runs) |
| **`Plasmid_analysis`** | *reads only its own inputs* |
| `LMM_all_results.csv` | LMM estimates per condition × replicon under each model |
| `replicon_VIF.csv` | Collinearity (VIF) between each replicon and the cargo genes |
| `replicon_cargo_cooccurrence.csv` | How often each replicon and cargo gene occur in the same isolate |
| `plasmid_best_hits_with_Kleborate.tsv` | ONT plasmid contigs with replicon hits and Kleborate cargo calls |
| `copy_number/copy_number_data.csv` | Plasmid copy number per isolate and replicon |
| `copy_number/copy_number_fitness_associations.csv` | Correlation between copy number and fitness per condition |
| `within_ST_matching/*.csv` | Within-ST propensity-score matching: matched fitness differences, diagnostics, covariate balance |

The revision scripts use 8 of the 25 `OriginalData/` files; the others are the full dataset behind the original analysis.

`Revision/Shared_input/` holds inputs used by more than one analysis:

| File | Contents | Used by |
|---|---|---|
| `Opacity_all_conditions_traits.csv` | Colony opacity per isolate and condition | `Opacity_analysis`, `Opacity_ST23` |
| `all_models_mean_sd_five_runs.csv` | Fitness-prediction model performance (mean ± SD over five runs) per condition | `Overfitting_analysis`, `Opacity_ST23` |

---

## Contact

This is a joint project developed by the **[Infectious Disease Epidemiology Lab](https://ide.kaust.edu.sa/) (KAUST)** and the **Banzhaf Lab (Newcastle University)**.

| Name | Role | Email |
|---|---|---|
| Danesh Moradigaravand | PI, KAUST | danesh.moradigaravand@kaust.edu.sa |
| Manuel Banzhaf | PI, Newcastle University | manuel.banzhaf@newcastle.ac.uk |
