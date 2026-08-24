# -*- coding: utf-8 -*-
"""
Article 3 — Pan-cancer ZP3 isoform PSI specificity comparison
=============================================
Compare PSI (percent spliced in) fingerprints of 7 ZP3 transcripts across TCGA pan-cancer, and test
Pan-cancer patterns of transcript PSI associations with immune features:

  1) PSI fingerprint heatmap (32 cancer types × 7 transcripts, median PSI)
  2) per-cancer transcript PSI × immune feature Spearman (32 × 7 × 7 = 1568 associations)
  3) Ecological analysis: cancer-level FL-PSI / RI-PSI median × ZP3-immune association strength (Avg_Rho)
  4) Isoform switch index log2(FL_PSI / RI_PSI) pan-cancer ranking

Data sources (all local):
  - zp3_isoform_proportions.csv (19131 samples × 7 transcript PSI)
  - TcgaTargetGtex_rsem_gene_tpm.gz (1.3G local TPM, immune gene extraction)
  - tcga_disease_map.json (sample barcode -> cancer type)
  - ensg_map.json (symbol -> Ensembl cache)

Output directory: zp3_psi_pancancer_results/
"""
import os
import sys
import json
import gzip
import time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))      # project root
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_pancancer_results")
os.makedirs(OUT, exist_ok=True)

DATA_TPM = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "phase1_knowledge_gap_filling", "TcgaTargetGtex_rsem_gene_tpm.gz")
PROP_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_isoform_proportions.csv")
DISEASE_MAP = os.path.join(ROOT, "output", "tcga_pancan", "tcga_disease_map.json")
ENSG_CACHE = os.path.join(ROOT, "output", "tcga_pancan", "ensg_map.json")
SUMMARY_CSV = os.path.join(ROOT, "output", "tcga_pancan", "tcga_pancan_cancer_summary.csv")

# Immune gene set identical to tcga_pancan_zp3_analysis.py
IMMUNE_GENE_SETS = {
    'M2_Macrophage': ['CD163', 'MSR1', 'MRC1', 'VSIG4', 'CD200R1', 'TGFB1', 'IL10',
                      'ARG1', 'MERTK', 'CLEC7A'],
    'T_cell_exhaustion': ['LAG3', 'TIGIT', 'HAVCR2', 'PDCD1', 'CTLA4', 'CD274',
                          'PDCD1LG2', 'BTLA', 'VSIR', 'IDO1', 'IDO2'],
    'Cytolytic_activity': ['GZMA', 'GZMB', 'PRF1', 'IFNG'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'ICOS', 'CD40LG'],
    'IFN_gamma': ['IFNG', 'STAT1', 'IRF1', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'CD274'],
    'Checkpoint': ['CD274', 'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2', 'BTLA', 'VSIR'],
    'Myeloid': ['CD68', 'CD163', 'CSF1R', 'ITGAM', 'CD14', 'LYZ', 'S100A8', 'S100A9'],
}
# Transcript semantic labels
TX_LABEL = {
    "ENST00000336517.8": "FL canonical (9-exon)",
    "ENST00000466960.5": "Retained-intron",
    "ENST00000394860.3": "5-exon truncated",
    "ENST00000394857.7": "Alt-exon isoform",
    "ENST00000467555.1": "Short isoform",
    "ENST00000416245.5": "Mid isoform",
    "ENST00000479793.5": "Alt-terminal isoform",
}
FL = "ENST00000336517.8"   # canonical full-length
RI = "ENST00000466960.5"   # retained-intron


def get_ensg_map(symbols):
    with open(ENSG_CACHE) as f:
        cache = json.load(f)
    return {s: cache.get(s) for s in symbols}


def read_target_genes(path, target_strip):
    """Stream-read target genes in TPM (target_strip: ensg without version -> symbol or ensg)."""
    rows = {}
    with gzip.open(path, "rt") as f:
        first = f.readline()
        samples = first.rstrip("\n").split("\t")[1:]
        print(f"    Total {len(samples)} samples")
        n = 0
        while True:
            lines = f.readlines(65536)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip("\n").split("\t")
                gid = parts[0].split(".")[0]
                if gid in target_strip:
                    rows[target_strip[gid]] = [float(x) for x in parts[1:]]
            n += len(lines)
            if n % 200000 == 0:
                print(f"    Scanned {n} gene rows...")
    return pd.DataFrame.from_dict(rows, orient="index", columns=samples)


def zscore_consensus_score(mat_gx_s, ensgs):
    """z-score consensus: take the mean after per-gene cross-sample standardization (consistent with the pan-cancer script). Returns sample-level Series."""
    sub = mat_gx_s.loc[ensgs]
    gene_mean = sub.mean(axis=1)
    gene_std = sub.std(axis=1)
    valid = gene_std > 0
    if not valid.any():
        return None
    z = (sub.loc[valid] - gene_mean[valid].values[:, None]) / gene_std[valid].values[:, None]
    return z.mean(axis=0)


def main():
    print("=== Article 3: Pan-cancer ZP3 isoform PSI specificity comparison ===\n")

    # ---- 1. PSI matrix ----
    print("1. Reading PSI proportion matrix ...")
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    print(f"   PSI matrix {psi.shape[0]} samples × {psi.shape[1]} transcripts")

    # ---- 2. Immune gene TPM (streaming) ----
    all_symbols = sorted({g for s in IMMUNE_GENE_SETS.values() for g in s})
    sym2ensg = get_ensg_map(all_symbols)
    unresolved = [s for s, e in sym2ensg.items() if not e]
    if unresolved:
        print(f"  !! Unresolved immune genes: {unresolved}")
    target_ids = [e for e in sym2ensg.values() if e]
    print(f"2. Streaming read real TPM ({len(target_ids)} immune genes)...")
    t0 = time.time()
    mat = read_target_genes(DATA_TPM, {e: e for e in target_ids})
    print(f"   Read complete {mat.shape[0]} genes × {mat.shape[1]} samples, elapsed {time.time()-t0:.1f}s")

    # ---- 3. Immune score (all TCGA tumor samples) ----
    samples = list(mat.columns)
    tcga_mask = [s.startswith("TCGA-") and s.split("-")[3].startswith("01") for s in samples]
    tcga_samples = [s for s, m in zip(samples, tcga_mask) if m]
    mat_t = mat[tcga_samples]
    print(f"3. TCGA tumor samples {len(tcga_samples)}, calculating 7-feature z-score consensus score ...")

    with open(DISEASE_MAP) as f:
        p2cancer = json.load(f)
    cancer_of = {}
    for s in tcga_samples:
        cancer_of[s] = p2cancer.get("-".join(s.split("-")[:3]), "UNKNOWN")

    score_mat = {}   # sample -> {feature: score}
    for set_name, syms in IMMUNE_GENE_SETS.items():
        ensgs = [sym2ensg[s] for s in syms if sym2ensg.get(s) in mat_t.index]
        sc = zscore_consensus_score(mat_t, ensgs)
        if sc is None:
            continue
        for s in tcga_samples:
            score_mat.setdefault(s, {})[set_name] = float(sc[s])
    print(f"   Scoring complete, {len(score_mat)} samples have scores")

    # ---- 4. join PSI + immune score + cancer type ----
    rec = []
    for s in tcga_samples:
        if s not in psi.index or s not in score_mat:
            continue
        c = cancer_of[s]
        if c == "UNKNOWN":
            continue
        r = {"Sample": s, "Cancer": c}
        for tx in psi.columns:
            r[tx] = float(psi.loc[s, tx])
        for fname, v in score_mat[s].items():
            r[f"score_{fname}"] = v
        rec.append(r)
    df = pd.DataFrame(rec)
    print(f"4. After merging {len(df)} samples (PSI × immune score × cancer type)")
    if len(df) == 0:
        print("!! Merge is empty, exiting"); sys.exit(1)
    df.to_csv(os.path.join(OUT, "psi_immune_joined_samples.csv"), index=False)

    # ---- 5. per-cancer transcript PSI × immune feature association ----
    print("5. per-cancer transcript PSI × immune feature Spearman ...")
    rows = []
    for c, grp in df.groupby("Cancer"):
        if len(grp) < 30:
            continue
        for tx in psi.columns:
            pv = grp[tx].values.astype(float)
            for fname in IMMUNE_GENE_SETS:
                sv = grp[f"score_{fname}"].values.astype(float)
                m = np.isfinite(pv) & np.isfinite(sv)
                if m.sum() < 20:
                    continue
                rho, p = stats.spearmanr(pv[m], sv[m])
                rows.append({
                    "Cancer": c, "Transcript": tx, "Tx_Label": TX_LABEL.get(tx, tx),
                    "Feature": fname, "Rho": round(float(rho), 4),
                    "P_value": float(p), "N": int(m.sum()),
                })
    corr = pd.DataFrame(rows)
    corr["FDR"] = 0.0
    for (c, tx), g in corr.groupby(["Cancer", "Transcript"]):
        if len(g) > 1:
            corr.loc[g.index, "FDR"] = _bh(g["P_value"].values)
        else:
            corr.loc[g.index, "FDR"] = g["P_value"].values
    corr.to_csv(os.path.join(OUT, "psi_immune_pancancer_correlation.csv"), index=False)
    print(f"   {len(corr)} associations, covering {corr['Cancer'].nunique()} cancer types")
    n_sig = (corr["FDR"] < 0.05).sum()
    print(f"   Significant (FDR<0.05): {n_sig} records")

    # FL transcript cross-cancer direction consistency
    fl = corr[corr["Transcript"] == FL]
    fl_pos_sig = ((fl["FDR"] < 0.05) & (fl["Rho"] > 0)).sum()
    fl_neg_sig = ((fl["FDR"] < 0.05) & (fl["Rho"] < 0)).sum()
    fl_all = len(fl)
    print(f"   FL associations: {fl_all} records, positive significant {fl_pos_sig} / negative significant {fl_neg_sig}")
    ri = corr[corr["Transcript"] == RI]
    ri_pos_sig = ((ri["FDR"] < 0.05) & (ri["Rho"] > 0)).sum()
    ri_neg_sig = ((ri["FDR"] < 0.05) & (ri["Rho"] < 0)).sum()
    print(f"   RI associations: {len(ri)} records, positive significant {ri_pos_sig} / negative significant {ri_neg_sig}")

    # ---- 6. cancer-level PSI fingerprint + ecology ----
    print("6. cancer-level PSI fingerprint and ecology analysis ...")
    fp = df.groupby("Cancer")[list(psi.columns)].median()
    fp["N"] = df.groupby("Cancer").size()
    fp = fp.reset_index().rename(columns={"index": "Cancer"})
    fp = fp.sort_values(FL, ascending=False)
    fp.to_csv(os.path.join(OUT, "psi_pancancer_fingerprint.csv"), index=False)

    # Switch index per-sample then take cancer-type median
    df["switch_index"] = np.log2((df[FL].clip(lower=1e-6)) / (df[RI].clip(lower=1e-6)))
    sw = df.groupby("Cancer")["switch_index"].median().reset_index()

    # Ecology: merge with pan-cancer Avg_Rho
    summ = pd.read_csv(SUMMARY_CSV)
    eco = fp.merge(summ, left_on="Cancer", right_on="Cancer_Code", how="inner")
    eco = eco.merge(sw, on="Cancer")
    eco.to_csv(os.path.join(OUT, "psi_pancancer_ecological.csv"), index=False)
    for tx, lab in [(FL, "FL canonical"), (RI, "Retained-intron")]:
        x = eco[tx].values.astype(float)
        y = eco["Avg_Rho"].values.astype(float)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() >= 5:
            rho, p = stats.spearmanr(x[m], y[m])
            print(f"   Ecological correlation: {lab}-PSI median × ZP3-immune Avg_Rho: ρ={rho:+.3f}, p={p:.3f}, n={m.sum()}")
        x2 = eco["switch_index"].values
        rho2, p2 = stats.spearmanr(x2[m], y[m])
        print(f"   Ecological correlation: switch_index × Avg_Rho: ρ={rho2:+.3f}, p={p2:.3f}")

    # ---- 7. Figure ----
    print("7. Plotting ...")
    plot_figure(fp, eco, corr)
    print("\n=== Done, outputs in", OUT, "===")


def _bh(pvals):
    """Benjamini-Hochberg FDR。"""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty(n)
    out[order] = q
    return out


def plot_figure(fp, eco, corr):
    fig = plt.figure(figsize=(15, 11))

    # (a) PSI fingerprint heatmap
    ax = fig.add_subplot(2, 2, 1)
    tx_cols = list(TX_LABEL.keys())
    data = fp.set_index("Cancer")[tx_cols]
    labels = [TX_LABEL[t] for t in tx_cols]
    sns.heatmap(data, cmap="YlOrRd", annot=False, fmt=".2f",
                linewidths=0.4, linecolor="white", cbar_kws={"label": "median PSI"},
                ax=ax, vmin=0, vmax=0.9)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
    ax.set_title(f"(a) ZP3 PSI fingerprint (n={len(fp)} cancers)", fontsize=12)

    # (b) Ecological: FL-PSI × Avg_Rho
    ax = fig.add_subplot(2, 2, 2)
    x = eco[FL].values.astype(float); y = eco["Avg_Rho"].values.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    rho, p = stats.spearmanr(x[m], y[m])
    ax.scatter(x[m], y[m], s=55, alpha=0.8, c="#378ADD", edgecolor="white")
    # Highlight glioma
    for cname, mk in [("GBM", "o"), ("LGG", "s")]:
        sub = eco[eco["Cancer"] == cname]
        if len(sub):
            ax.scatter(sub[FL], sub["Avg_Rho"], s=110, c="#C00000", marker=mk,
                       edgecolor="black", zorder=5, label=cname)
    xx = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
    if m.sum() >= 5:
        slope, intercept = np.polyfit(x[m], y[m], 1)
        ax.plot(xx, slope * xx + intercept, "k--", lw=1, alpha=0.6)
    ax.set_xlabel("Median PSI of FL canonical transcript")
    ax.set_ylabel("ZP3-immune Avg_Rho (pancancer)")
    ax.set_title(f"(b) Ecological: FL-PSI × immune assoc\nρ={rho:+.3f}, p={p:.3f}, n={m.sum()}",
                 fontsize=11)
    ax.legend(fontsize=8)
    ax.axhline(0, color="grey", lw=0.6, ls=":")

    # (c) Ecological: RI-PSI × Avg_Rho
    ax = fig.add_subplot(2, 2, 3)
    x = eco[RI].values.astype(float); y = eco["Avg_Rho"].values.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    rho, p = stats.spearmanr(x[m], y[m])
    ax.scatter(x[m], y[m], s=55, alpha=0.8, c="#1D9E75", edgecolor="white")
    for cname, mk in [("GBM", "o"), ("LGG", "s")]:
        sub = eco[eco["Cancer"] == cname]
        if len(sub):
            ax.scatter(sub[RI], sub["Avg_Rho"], s=110, c="#C00000", marker=mk,
                       edgecolor="black", zorder=5, label=cname)
    if m.sum() >= 5:
        xx = np.linspace(np.nanmin(x[m]), np.nanmax(x[m]), 50)
        slope, intercept = np.polyfit(x[m], y[m], 1)
        ax.plot(xx, slope * xx + intercept, "k--", lw=1, alpha=0.6)
    ax.set_xlabel("Median PSI of retained-intron transcript")
    ax.set_ylabel("ZP3-immune Avg_Rho (pancancer)")
    ax.set_title(f"(c) Ecological: RI-PSI × immune assoc\nρ={rho:+.3f}, p={p:.3f}, n={m.sum()}",
                 fontsize=11)
    ax.legend(fontsize=8)
    ax.axhline(0, color="grey", lw=0.6, ls=":")

    # (d) FL transcript PSI × immune features per-cancer correlation heatmap
    ax = fig.add_subplot(2, 2, 4)
    fl = corr[corr["Transcript"] == FL]
    piv = fl.pivot_table(index="Cancer", columns="Feature", values="Rho")
    piv = piv.reindex(fp["Cancer"])
    sns.heatmap(piv, cmap="RdBu_r", center=0, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "ρ (FL-PSI × feature)"}, ax=ax,
                vmin=-0.4, vmax=0.4)
    ax.set_title("(d) FL-PSI × immune feature, per cancer", fontsize=12)
    ax.tick_params(axis="y", labelsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_zp3_psi_pancancer.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
