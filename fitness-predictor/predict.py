"""
predict.py — Predict fitness for a new strain from its genome FASTA.

Usage:
    fitness-predict -c <condition> -p <input.fasta>
    fitness-predict --download-data <directory>

Examples:
    fitness-predict --download-data ~/fitness_data
    fitness-predict -c LB -p query.fasta
    fitness-predict -c Citricacid -p query.fasta --identity 90 --coverage 70
    fitness-predict -c Ampicillin128ugml -p query.fasta --explain
    fitness-predict -c Ampicillin128ugml -p query.fasta --explain --top-n 30

Options:
    -c, --condition       Condition name (must match a trained model in ./models/)
    -p, --fasta           Input genome FASTA file
    --identity            Min BLAST % identity   [default: 95]
    --coverage            Min query coverage %   [default: 80]
    --explain             Compute SHAP values and report top contributing genes
    --top-n N             Number of top genes to report with --explain [default: 20]
    --list                List available trained models and exit
    --download-data DIR   Download models and reference to DIR and exit
"""

import argparse
import os
import sys
import subprocess
import tempfile
import numpy as np
import pandas as pd
import joblib


DATA_DIR   = os.environ.get("FITNESS_DATA_DIR", os.getcwd())
MODEL_DIR  = os.path.join(DATA_DIR, "models")
PAN_DIR    = os.path.join(DATA_DIR, "input_pan_genome")
REF_FA     = os.path.join(DATA_DIR, "pan_genome_reference.fa")
ROARY_CSV  = os.path.join(DATA_DIR, "gene_annotations.csv")


# ──────────────────────────────────────────────────────────────
# Gene annotations
# ──────────────────────────────────────────────────────────────

def load_annotations():
    if not os.path.exists(ROARY_CSV):
        return {}
    df = pd.read_csv(ROARY_CSV, usecols=["Gene", "Annotation"])
    return dict(zip(df["Gene"], df["Annotation"]))


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
# SHAP
# ──────────────────────────────────────────────────────────────

def compute_shap(model, feature_cols, x_arr):
    """
    Compute per-gene SHAP values for a single sample using the NGBoost
    TreeExplainer approach.  Background = all-zeros (no genes present).

    Returns base_value (fitness predicted with no genes) and a dict
    {gene: shap_value} for genes present in the sample.
    """
    try:
        import shap
    except ImportError:
        sys.exit("[ERROR] SHAP not installed. Run: pip install shap")

    n  = len(feature_cols)
    X  = x_arr.reshape(1, -1)
    bg = np.zeros((1, n))

    base_value = float(model.predict(bg)[0])

    bm       = model.base_models
    lr       = model.learning_rate
    scalings = model.scalings
    col_idxs = model.col_idxs

    sv_total = np.zeros(n)
    for step_trees, sc, col_idx in zip(bm, scalings, col_idxs):
        te      = shap.TreeExplainer(step_trees[0], data=bg[:, col_idx])
        sv_step = np.zeros(n)
        sv_step[col_idx] = te.shap_values(X[:, col_idx])[0]
        sv_total -= sc * lr * sv_step   # NGBoost subtracts gradients

    present   = np.where(x_arr == 1.0)[0]
    gene_shap = {feature_cols[i]: round(float(sv_total[i]), 6) for i in present}
    return base_value, gene_shap


# ──────────────────────────────────────────────────────────────
# Prediction
# ──────────────────────────────────────────────────────────────

def predict(condition, fasta, identity=95.0, coverage=80.0, explain=False, top_n=20):
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

    result = {
        "condition":    condition,
        "fasta":        fasta,
        "genes_hit":    int(x_vec.sum()),
        "total_genes":  len(feature_cols),
        "fitness_mean": round(float(y_mean), 5),
        "ci95_lower":   round(float(lower),  5),
        "ci95_upper":   round(float(upper),  5),
    }

    if explain:
        print("Computing SHAP values…")
        x_arr       = x_vec.values
        annotations = load_annotations()
        base_value, gene_shap = compute_shap(model, list(feature_cols), x_arr)
        sorted_genes = sorted(gene_shap.items(), key=lambda kv: abs(kv[1]), reverse=True)
        result["shap_base_value"] = round(base_value, 5)
        result["shap_top_genes"]  = [
            {"gene": g, "shap": v, "annotation": annotations.get(g, "")}
            for g, v in sorted_genes[:top_n]
        ]

    return result


# ──────────────────────────────────────────────────────────────
# Data download
# ──────────────────────────────────────────────────────────────

RELEASE_BASE = "https://zenodo.org/records/19819400/files"

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

    # Download gene annotations (for --explain)
    roary_url  = f"{RELEASE_BASE}/gene_annotations.csv"
    roary_dest = os.path.join(dest_dir, "gene_annotations.csv")
    if not os.path.exists(roary_dest):
        print(f"Downloading gene_annotations.csv…")
        urllib.request.urlretrieve(roary_url, roary_dest)
    else:
        print(f"gene_annotations.csv already exists, skipping.")

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
    parser.add_argument("--explain", action="store_true",
                        help="Compute SHAP values and report top contributing genes")
    parser.add_argument("--top-n", type=int, default=20, metavar="N",
                        help="Number of top genes to show with --explain [default: 20]")
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
        explain   = args.explain,
        top_n     = args.top_n,
    )

    print()
    print(f"Condition      : {result['condition']}")
    print(f"FASTA          : {result['fasta']}")
    print(f"Genes detected : {result['genes_hit']} / {result['total_genes']}")
    print(f"Fitness (mean) : {result['fitness_mean']:.4f}")
    print(f"95% CI         : [{result['ci95_lower']:.4f}, {result['ci95_upper']:.4f}]")

    if args.explain and "shap_top_genes" in result:
        print()
        print(f"SHAP (baseline fitness with no genes: {result['shap_base_value']:.4f})")
        print(f"  {'Gene':<28} {'SHAP':>10}  {'Effect':<10}  Annotation")
        print(f"  {'-'*28} {'-'*10}  {'-'*10}  ----------")
        for entry in result["shap_top_genes"]:
            v      = entry["shap"]
            effect = "promotes" if v > 0 else "reduces" if v < 0 else "neutral"
            annot  = entry.get("annotation", "")[:60]
            print(f"  {entry['gene']:<28} {v:>+10.5f}  {effect:<10}  {annot}")

    print()

    base_name = os.path.splitext(args.fasta)[0]
    if base_name.endswith(".fasta") or base_name.endswith(".fa") or base_name.endswith(".fna"):
        base_name = os.path.splitext(base_name)[0]

    # Prediction summary CSV
    pred_csv = f"{base_name}_{args.condition}_prediction.csv"
    pd.DataFrame([{
        "condition":    result["condition"],
        "fasta":        result["fasta"],
        "genes_hit":    result["genes_hit"],
        "total_genes":  result["total_genes"],
        "fitness_mean": result["fitness_mean"],
        "ci95_lower":   result["ci95_lower"],
        "ci95_upper":   result["ci95_upper"],
    }]).to_csv(pred_csv, index=False)
    print(f"Result saved   : {pred_csv}")

    # SHAP CSV
    if "shap_top_genes" in result:
        shap_csv = f"{base_name}_{args.condition}_shap.csv"
        rows = []
        for entry in result["shap_top_genes"]:
            v = entry["shap"]
            rows.append({
                "gene":       entry["gene"],
                "shap":       v,
                "effect":     "promotes" if v > 0 else "reduces" if v < 0 else "neutral",
                "annotation": entry.get("annotation", ""),
            })
        pd.DataFrame(rows).to_csv(shap_csv, index=False)
        print(f"SHAP saved     : {shap_csv}")


if __name__ == "__main__":
    main()
