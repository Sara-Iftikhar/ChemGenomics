# Bacterial Fitness Predictor

Predict the absolute fitness of a bacterial strain under different growth conditions from a genome FASTA file using NGBoost probabilistic regression.

## Installation

### Option 1 — Conda

Requires [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

```bash
conda install -c gzhoubioinf -c bioconda -c conda-forge fitness-predictor
```

Then download models and reference data (one time):

```bash
fitness-predict --download-data ~/fitness_data
cd ~/fitness_data
```

### Option 2 — Docker

Requires [Docker](https://docs.docker.com/get-docker/). Models and reference are bundled in the image — no data download needed.

```bash
docker pull gzhoubioinf09/fitness-predictor:v0.1.0

# List available conditions
docker run --rm gzhoubioinf09/fitness-predictor:v0.1.0 --list

# Predict fitness (mount your directory to access FASTA and save output)
docker run --rm -v $(pwd):/data gzhoubioinf09/fitness-predictor:v0.1.0 \
    -c <condition> -p /data/genome.fasta
```

## Usage

> For conda users: run commands from `~/fitness_data` (where models are downloaded).

```bash
# List all available conditions (160 total)
fitness-predict --list

# Predict fitness for a strain (FASTA can be any path)
fitness-predict -c <condition> -p /path/to/genome.fasta

# Predict with gzipped FASTA
fitness-predict -c <condition> -p /path/to/genome.fasta.gz

# Custom BLAST thresholds (default: identity=95, coverage=80)
fitness-predict -c <condition> -p /path/to/genome.fasta --identity 80 --coverage 60
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

- **Fitness > 1** — strain grows better under this condition than in the control (LB)
- **Fitness ≈ 1** — strain grows similarly to the control
- **Fitness < 1** — strain grows worse than in the control (LB)

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
fitness-predictor/
├── predict.py              # CLI prediction script
├── environment.yml         # Conda environment
├── requirements.txt        # Pip dependencies
├── pan_genome_reference.fa # BLAST reference (pan-genome genes)
└── models/                 # 160 trained NGBoost models (.joblib)
```
