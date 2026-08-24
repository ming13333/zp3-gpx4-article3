#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ZP3 real isoform quantification (based on TCGA TARGET GTEx isoform TPM)
=========================================================
Data: local TcgaTargetGtex_rsem_isoform_tpm.gz (4.3 GB)
Method:
  1. Stream-read, keeping only ZP3 candidate transcripts (8 Ensembl plus-strand transcripts).
  2. Convert log2(TPM) -> TPM (relative quantity), and calculate the proportion of each transcript to ZP3 total expression per sample.
  3. Use GDC mapping (tcga_disease_map.json) for TCGA cancer-type grouping;
     Use GTEX_phenotype.gz to distinguish GTEx normal samples.
  4. Compare isoform proportion differences in tumor vs normal (Mann-Whitney U + BH FDR).

Output:
  - zp3_isoform_proportions.csv: per-sample proportions of each transcript
  - zp3_isoform_tumor_vs_normal.csv: tumor vs normal differential test
  - zp3_isoform_by_cancer.csv: median isoform proportions by cancer type
  - fig_zp3_isoform_proportions.png: visualization
"""
import os, sys, json, gzip, time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(BASE)))
DATA_ISO = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling",
                        "TcgaTargetGtex_rsem_isoform_tpm.gz")
DATA_PHENO = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling",
                           "GTEX_phenotype.gz")
GDC_MAP = os.path.join(ROOT, "output", "tcga_pancan", "tcga_disease_map.json")
OUT_DIR = os.path.join(ROOT, "article3", "results")

# ZP3 plus-strand candidate transcripts (Ensembl GRCh38, obtained from overlap query)
ZP3_TRANSCRIPTS = {
    'ENST00000336517', 'ENST00001135277', 'ENST00000394857',
    'ENST00000416245', 'ENST00000394860', 'ENST00000466960',
    'ENST00000479793', 'ENST00000467555'
}

# ---------------------------------------------------------------------------
# 1. Stream-read isoform TPM (only ZP3 candidate transcripts)
# ---------------------------------------------------------------------------
def read_zp3_isoforms(path, target_prefixes, chunk=5000):
    """Stream-read the gz isoform matrix, keeping only target_prefixes (without version number).
    Returns DataFrame: index=transcript ID (with version number), columns=sample id."""
    rows = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        n = 0
        while True:
            lines = f.readlines(chunk)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip("\n").split("\t")
                tid = parts[0]
                base = tid.split(".")[0]
                if base in target_prefixes:
                    rows[tid] = [float(x) for x in parts[1:]]
            n += len(lines)
            if n % 20000 == 0:
                print(f"    Scanned {n} transcript lines...")
    df = pd.DataFrame.from_dict(rows, orient="index", columns=samples)
    return df


# ---------------------------------------------------------------------------
# 2. Load GDC cancer type mapping + GTEx phenotypes
# ---------------------------------------------------------------------------
def load_gdc_map(path):
    with open(path) as f:
        return json.load(f)


def load_gtex_pheno(path):
    """GTEX_phenotype.gz is tsv, containing Sample column and _primary_site etc."""
    df = pd.read_csv(path, sep="\t", compression="gzip")
    # Sample column name may be 'Sample' or the first column
    sample_col = 'Sample' if 'Sample' in df.columns else df.columns[0]
    df = df.rename(columns={sample_col: 'sample_id'})
    return df


# ---------------------------------------------------------------------------
# 3. Main workflow
# ---------------------------------------------------------------------------
def main():
    print("=== ZP3 real isoform quantification (based on TCGA TARGET GTEx isoform TPM) ===\n")
    if not os.path.exists(DATA_ISO):
        print(f"!! Isoform data not found: {DATA_ISO}"); sys.exit(1)

    # 3.1 Read isoforms
    print("1. Stream-read isoform TPM (ZP3 candidate transcripts only)...")
    t0 = time.time()
    iso = read_zp3_isoforms(DATA_ISO, ZP3_TRANSCRIPTS)
    print(f"   Found {iso.shape[0]} ZP3 transcripts × {iso.shape[1]} samples,"
          f"took {time.time()-t0:.1f}s")
    if iso.empty:
        print("!! No ZP3 transcripts found, exiting"); sys.exit(1)

    # 3.2 Convert log2(TPM) -> TPM (relative quantity)
    print("\n2. Convert log2(TPM) -> TPM and calculate per-sample isoform proportions...")
    tpm = 2 ** iso  # Inverse transformation of log2(TPM), yielding relative TPM
    # Avoid negative values (theoretically 2^log2(TPM) should be positive, but the minimum value -9.9658 yields extremely small values)
    tpm = tpm.clip(lower=0)
    # Total ZP3 expression per sample
    total_zp3 = tpm.sum(axis=0)
    # Proportion of each transcript per sample (set proportions to NaN for samples with total ZP3 = 0).
    # Note: do not use prop.where(cond) -- pandas aligns the Series condition with the row index
    # rather than column index, causing all NaN; changed to explicit per-column assignment.
    prop = tpm.div(total_zp3, axis=1)
    prop.loc[:, total_zp3 <= 0] = np.nan
    print(f"   Number of samples with total ZP3 TPM=0: {(total_zp3 == 0).sum()}")

    # Save proportion matrix
    prop_path = os.path.join(OUT_DIR, "zp3_isoform_proportions.csv")
    prop.T.to_csv(prop_path)
    print(f"   Proportion matrix saved: {prop_path}")

    # 3.3 Sample classification
    samples = list(prop.columns)
    # TCGA tumor: starts with TCGA- and sample segment starts with 01
    tcga_tumor_mask = [s.startswith("TCGA-") and s.split("-")[3].startswith("01")
                       for s in samples]
    tcga_tumor_samples = [s for s, m in zip(samples, tcga_tumor_mask) if m]
    # GTEx normal
    gtex_samples = [s for s in samples if s.startswith("GTEX-")]
    print(f"\n3. Sample classification:")
    print(f"   TCGA tumor samples: {len(tcga_tumor_samples)}")
    print(f"   GTEx normal samples: {len(gtex_samples)}")

    # 3.4 Tumor vs normal: isoform proportion differences
    print("\n4. Tumor vs normal (GTEx) isoform proportion differences...")
    tumor_prop = prop[tcga_tumor_samples].T
    normal_prop = prop[gtex_samples].T
    results = []
    for tid in prop.index:
        t_vals = tumor_prop[tid].dropna()
        n_vals = normal_prop[tid].dropna()
        if len(t_vals) < 10 or len(n_vals) < 10:
            continue
        u, p = stats.mannwhitneyu(t_vals, n_vals, alternative='two-sided')
        r = 1 - (2 * u) / (len(t_vals) * len(n_vals))  # rank-biserial
        results.append({
            'Transcript': tid,
            'Tumor_median': t_vals.median(),
            'Normal_median': n_vals.median(),
            'Tumor_mean': t_vals.mean(),
            'Normal_mean': n_vals.mean(),
            'MannWhitney_p': p,
            'Effect_r': r,
            'Tumor_n': len(t_vals),
            'Normal_n': len(n_vals),
        })
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(res_df['MannWhitney_p'], method='fdr_bh')
        res_df['FDR'] = fdr
        res_df = res_df.sort_values('FDR')
    res_path = os.path.join(OUT_DIR, "zp3_isoform_tumor_vs_normal.csv")
    res_df.to_csv(res_path, index=False)
    print(f"   results saved: {res_path}")
    print(res_df.to_string(index=False))

    # 3.5 By cancer type analysis (TCGA tumor)
    print("\n5. Analyzing isoform proportions by cancer type...")
    if os.path.exists(GDC_MAP):
        p2cancer = load_gdc_map(GDC_MAP)
        participant_of = {"-".join(s.split("-")[:3]) for s in tcga_tumor_samples}
        cancer_of = {s: p2cancer.get("-".join(s.split("-")[:3]), "UNKNOWN")
                     for s in tcga_tumor_samples}
        cancers = {}
        for s in tcga_tumor_samples:
            c = cancer_of[s]
            if c != "UNKNOWN":
                cancers.setdefault(c, []).append(s)
        print(f"   covers {len(cancers)} cancer types")

        cancer_records = []
        for cancer, sams in cancers.items():
            if len(sams) < 20:
                continue
            sub = prop[sams].T
            for tid in prop.index:
                vals = sub[tid].dropna()
                if len(vals) < 10:
                    continue
                cancer_records.append({
                    'Cancer': cancer, 'Transcript': tid,
                    'Median_prop': vals.median(), 'Mean_prop': vals.mean(),
                    'N': len(vals),
                })
        cancer_df = pd.DataFrame(cancer_records)
        cancer_path = os.path.join(OUT_DIR, "zp3_isoform_by_cancer.csv")
        cancer_df.to_csv(cancer_path, index=False)
        print(f"   results saved: {cancer_path}")
    else:
        print(f"   (skipped: GDC mapping cache not found: {GDC_MAP})")

    # 3.6 Visualization
    print("\n6. Generating visualizations...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # (a) Tumor vs Normal box plot
    plot_data = []
    for tid in prop.index:
        t_vals = tumor_prop[tid].dropna()
        n_vals = normal_prop[tid].dropna()
        for v in t_vals:
            plot_data.append({'Transcript': tid.split('.')[0], 'Group': 'Tumor', 'Proportion': v})
        for v in n_vals:
            plot_data.append({'Transcript': tid.split('.')[0], 'Group': 'Normal', 'Proportion': v})
    plot_df = pd.DataFrame(plot_data)
    if not plot_df.empty:
        sns.boxplot(data=plot_df, x='Transcript', y='Proportion', hue='Group', ax=axes[0])
        axes[0].set_title('ZP3 Isoform Proportions: Tumor vs Normal')
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].set_ylabel('Proportion of total ZP3')

    # (b) Transcript expression (log2 TPM) violin plot
    log_data = []
    for tid in iso.index:
        for s in tcga_tumor_samples[:500]:  # Sample to avoid overcrowding the plot
            v = iso.loc[tid, s]
            if np.isfinite(v):
                log_data.append({'Transcript': tid.split('.')[0], 'log2TPM': v})
    log_df = pd.DataFrame(log_data)
    if not log_df.empty:
        sns.violinplot(data=log_df, x='Transcript', y='log2TPM', ax=axes[1])
        axes[1].set_title('ZP3 Isoform Expression Distribution (TCGA Tumor)')
        axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    fig_path = os.path.join(OUT_DIR, "fig_zp3_isoform_proportions.png")
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"   Figure saved: {fig_path}")

    print("\n=== analysis complete ===")


if __name__ == "__main__":
    main()
