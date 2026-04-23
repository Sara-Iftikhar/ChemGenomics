# Bacterial Fitness Predictor

Predict the absolute fitness of a bacterial strain under different growth conditions from a genome FASTA file using NGBoost probabilistic regression.

## Requirements

- [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Installation

```bash
conda install -c gzhoubioinf -c bioconda -c conda-forge fitness-predictor
```

## Setup (one time)

```bash
# Download models and reference data
fitness-predict --download-data ~/fitness_data

# Go to the data directory
cd ~/fitness_data
```

## Usage

```bash
cd ~/fitness_data

# Predict fitness for a strain (FASTA can be any path)
fitness-predict -c Doripenem0125ugml -p /path/to/genome.fasta

# Predict with gzipped FASTA
fitness-predict -c LB -p /path/to/genome.fasta.gz

# List all available conditions (160 total)
fitness-predict --list

# Custom BLAST thresholds (default: identity=95, coverage=80)
fitness-predict -c Citricacid -p /path/to/genome.fasta --identity 80 --coverage 60
```

### Output

```
Condition      : Doripenem0125ugml
FASTA          : genome.fasta
Genes detected : 412 / 4751
Fitness (mean) : 1.0842
95% CI         : [0.8901, 1.2783]

Result saved   : genome_Doripenem0125ugml_prediction.json
```

- **Fitness > 1** — strain grows better than average under this condition
- **Fitness ≈ 1** — strain grows similarly to average
- **Fitness < 1** — strain grows worse than average

## Accepted input formats

| Extension | Supported |
|---|---|
| `.fasta`, `.fa`, `.fna` | Yes |
| `.fasta.gz`, `.fa.gz`, `.fna.gz` | Yes (decompressed automatically) |
| `.fasta.bz2`, `.fasta.zst` | No |

## How it works

1. The input genome is BLASTed against a pan-genome reference (`pan_genome_reference.fa`)
2. Gene presence/absence is determined from BLAST hits (identity ≥ 95%, coverage ≥ 80% by default)
3. A pre-trained NGBoost model predicts fitness and a 95% confidence interval from the gene vector

Models were trained on 160 conditions using GWAS-significant genes (Benjamini-Hochberg p < 0.05) as features and Spearman rho as the validation metric (80/20 train/test split).

## File structure

```
fitness_predictor/
├── predict.py              # CLI prediction script
├── environment.yml         # Conda environment
├── requirements.txt        # Pip dependencies
├── pyproject.toml          # Installs fitness-predict command
├── pan_genome_reference.fa # BLAST reference (pan-genome genes)
└── models/                 # 160 trained NGBoost models (.joblib)
```
