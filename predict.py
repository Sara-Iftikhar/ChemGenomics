"""
predict.py — Predict fitness for a new strain from its genome FASTA.

Usage:
    fitness-predict -c <condition> -p <input.fasta>
    fitness-predict --download-data <directory>

Examples:
    fitness-predict --download-data ~/fitness_data
    fitness-predict -c LB -p query.fasta
    fitness-predict -c Citricacid -p query.fasta --identity 90 --coverage 70

Options:
    -c, --condition       Condition name (must match a trained model in ./models/)
    -p, --fasta           Input genome FASTA file
    --identity            Min BLAST % identity   [default: 95]
    --coverage            Min query coverage %   [default: 80]
    --list                List available trained models and exit
    --download-data DIR   Download models and reference to DIR and exit
"""

import argparse
import os
import sys
import subprocess
import tempfile
import json
import numpy as np
import pandas as pd
import joblib


DATA_DIR  = os.environ.get("FITNESS_DATA_DIR", os.getcwd())
MODEL_DIR = os.path.join(DATA_DIR, "models")
PAN_DIR   = os.path.join(DATA_DIR, "input_pan_genome")
REF_FA    = os.path.join(DATA_DIR, "pan_genome_reference.fa")


# ──────────────────────────────────────────────────────────────
# BLAST helpers
# ──────────────────────────────────────────────────────────────

def check_blast():
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        sys.exit(
            "[ERROR] BLAST not found. Install with:\n"
            "  conda install -c bioconda blast\n"
            "  or: brew install blast"
        )


def make_blast_db(ref_fa):
    db_dir  = os.path.join(DATA_DIR, ".blastdb")
    db_path = os.path.join(db_dir, "pan_ref")
    flag    = db_path + ".nhr"

    if not os.path.exists(flag):
        os.makedirs(db_dir, exist_ok=True)
        print("Building BLAST database (one-time)…")
        subprocess.run(
            ["makeblastdb", "-in", ref_fa, "-dbtype", "nucl", "-out", db_path],
            check=True, capture_output=True
        )
    return db_path


def run_blast(query_fa, db_path, identity, coverage):
    """
    Run blastn and return set of reference gene names that pass thresholds.
    Uses outfmt 6: qseqid sseqid pident length slen
    Handles gzipped FASTA automatically.
    """
    # Decompress if gzipped
    if query_fa.endswith(".gz"):
        import gzip, shutil
        tmp_fa = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
        with gzip.open(query_fa, 'rb') as f_in, open(tmp_fa.name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        query_fa = tmp_fa.name
        gz_tmp = tmp_fa.name
    else:
        gz_tmp = None

    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        out_file = tmp.name

    subprocess.run(
        [
            "blastn",
            "-query",    query_fa,
            "-db",       db_path,
            "-out",      out_file,
            "-outfmt",   "6 qseqid sseqid pident length slen",
            "-perc_identity", str(identity),
            "-evalue",   "1e-5",
            "-num_threads", "4",
        ],
        check=True, capture_output=True
    )

    hits = set()
    with open(out_file) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            _qid, sseqid, pident, length, slen = parts[:5]
            # Coverage = fraction of the reference GENE covered by the alignment
            cov = 100.0 * float(length) / float(slen)
            if float(pident) >= identity and cov >= coverage:
                hits.add(sseqid.strip())

    os.unlink(out_file)
    if gz_tmp:
        os.unlink(gz_tmp)
    return hits


# ──────────────────────────────────────────────────────────────
# Feature vector
# ──────────────────────────────────────────────────────────────

def build_feature_vector(hits, feature_cols):
    vec = pd.Series(0, index=feature_cols, dtype=float)
    for gene in hits:
        if gene in vec.index:
            vec[gene] = 1.0
    return vec


# ──────────────────────────────────────────────────────────────
# Prediction
# ──────────────────────────────────────────────────────────────

def predict(condition, fasta, identity=95.0, coverage=80.0):
    model_path = os.path.join(MODEL_DIR, f"{condition}_ngboost.joblib")

    if not os.path.exists(model_path):
        sys.exit(
            f"[ERROR] No trained model for '{condition}'.\n"
            f"  Run: fitness-predict --list to see available conditions.\n"
            f"  DATA_DIR: {DATA_DIR}"
        )

    bundle       = joblib.load(model_path)
    model        = bundle["model"]
    feature_cols = bundle["feature_cols"]

    # BLAST
    check_blast()
    db_path = make_blast_db(REF_FA)

    print(f"Running BLAST (identity≥{identity}%, coverage≥{coverage}%)…")
    hits  = run_blast(fasta, db_path, identity, coverage)
    x_vec = build_feature_vector(hits, feature_cols)
    print(f"  Genes detected: {int(x_vec.sum())} / {len(feature_cols)}")
    X     = pd.DataFrame([x_vec])

    y_mean = model.predict(X)[0]
    dist   = model.pred_dist(X)
    scale  = dist.scale[0]
    lower  = y_mean - 1.96 * scale
    upper  = y_mean + 1.96 * scale

    return {
        "condition":    condition,
        "fasta":        fasta,
        "genes_hit":    int(x_vec.sum()),
        "total_genes":  len(feature_cols),
        "fitness_mean": round(float(y_mean), 5),
        "ci95_lower":   round(float(lower),  5),
        "ci95_upper":   round(float(upper),  5),
    }


# ──────────────────────────────────────────────────────────────
# Data download
# ──────────────────────────────────────────────────────────────

RELEASE_BASE = "https://zenodo.org/records/19707430/files"

def download_data(dest_dir):
    import urllib.request, tarfile
    os.makedirs(dest_dir, exist_ok=True)

    # Download pan-genome reference
    ref_url  = f"{RELEASE_BASE}/pan_genome_reference.fa"
    ref_dest = os.path.join(dest_dir, "pan_genome_reference.fa")
    if not os.path.exists(ref_dest):
        print(f"Downloading pan_genome_reference.fa…")
        urllib.request.urlretrieve(ref_url, ref_dest)
    else:
        print(f"pan_genome_reference.fa already exists, skipping.")

    # Download models tarball
    models_url  = f"{RELEASE_BASE}/models.tar.gz"
    models_dest = os.path.join(dest_dir, "models.tar.gz")
    models_dir  = os.path.join(dest_dir, "models")
    if not os.path.isdir(models_dir):
        print(f"Downloading models.tar.gz…")
        urllib.request.urlretrieve(models_url, models_dest)
        print(f"Extracting models…")
        with tarfile.open(models_dest) as tar:
            tar.extractall(dest_dir)
        os.unlink(models_dest)
    else:
        print(f"models/ already exists, skipping.")

    print(f"\nData ready in: {dest_dir}")
    print(f"Run predictions from that directory:")
    print(f"  cd {dest_dir}")
    print(f"  fitness-predict -c <condition> -p /path/to/genome.fasta")


# ──────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────

def list_models():
    if not os.path.isdir(MODEL_DIR):
        print("No models trained yet. Run: python train.py -c <condition>")
        return
    models = [f.replace("_ngboost.joblib", "")
              for f in os.listdir(MODEL_DIR) if f.endswith("_ngboost.joblib")]
    if models:
        print("Trained models:")
        for m in sorted(models):
            print(f"  {m}")
    else:
        print("No trained models found.")


def main():
    parser = argparse.ArgumentParser(
        description="Predict fitness from a genome FASTA using a trained NGBoost model"
    )
    parser.add_argument("-c", "--condition", metavar="COND",
                        help="Condition name (e.g. LB, Citricacid)")
    parser.add_argument("-p", "--fasta", metavar="FILE",
                        help="Input genome FASTA file")
    parser.add_argument("--identity", type=float, default=95.0,
                        help="Min BLAST %% identity [default: 95]")
    parser.add_argument("--coverage", type=float, default=80.0,
                        help="Min query coverage %% [default: 80]")
    parser.add_argument("--list", action="store_true",
                        help="List available trained models and exit")
    parser.add_argument("--download-data", metavar="DIR",
                        help="Download models and reference data to DIR and exit")
    args = parser.parse_args()

    if args.download_data:
        download_data(args.download_data)
        return

    if args.list:
        list_models()
        return

    if not args.condition or not args.fasta:
        parser.error("Both -c/--condition and -p/--fasta are required")

    if not os.path.isfile(args.fasta):
        sys.exit(f"[ERROR] FASTA file not found: {args.fasta}")

    result = predict(
        condition = args.condition,
        fasta     = args.fasta,
        identity  = args.identity,
        coverage  = args.coverage,
    )

    print()
    print(f"Condition      : {result['condition']}")
    print(f"FASTA          : {result['fasta']}")
    print(f"Genes detected : {result['genes_hit']} / {result['total_genes']}")
    print(f"Fitness (mean) : {result['fitness_mean']:.4f}")
    print(f"95% CI         : [{result['ci95_lower']:.4f}, {result['ci95_upper']:.4f}]")
    print()

    # Also write JSON result
    out_json = os.path.splitext(args.fasta)[0] + f"_{args.condition}_prediction.json"
    with open(out_json, "w") as fh:
        import json
        json.dump(result, fh, indent=2)
    print(f"Result saved   : {out_json}")


if __name__ == "__main__":
    main()
