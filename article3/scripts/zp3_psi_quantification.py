#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Article 3: ZP3 transcript PSI quantification + immune association
============================================
Based on the existing Toil isoform TPM (TcgaTargetGtex_rsem_isoform_tpm.gz, 4.2G)
Perform PSI (Percent Spliced In) quantification for the 7 ZP3 transcripts:
  PSI_t = TPM(transcript t) / Σ_all ZP3 transcripts TPM

Analysis:
  1. All TCGA tumors vs GTEx normal: PSI differences (tumor_vs_normal.csv already exists, focusing here on GBM/LGG)
  2. GBM(153) vs LGG(509): PSI comparison between tumor types
  3. **Transcript-level immune association**: PSI of each transcript in GBM/LGG vs immune features
     (M2 / Treg / Checkpoint z-score consensus) Spearman correlation — to answer
     "which transcript drives the ZP3-immunosuppression association"
  4. Key isoform: ENST00000394860.3 (5-exon truncated, tumor-enriched 22-fold)
     PSI-specific association

Data:
  - zp3_isoform_proportions.csv (19131 samples × 7 transcripts, produced by real_quant)
  - h2_bulk/TCGA.GBM.sampleMap and TCGA.LGG.sampleMap (21 immune genes log2TPM)

Outputs:
  psi_status_by_transcript.csv  —— GBM/LGG per-transcript PSI median + GBM vs LGG MWU
  psi_immune_correlation.csv    —— Spearman correlation of each transcript PSI × immune feature
  fig_zp3_psi_immune.png        —— visualization
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
H2 = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "h2_bulk")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_results")
os.makedirs(OUT, exist_ok=True)

PROP = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_isoform_proportions.csv")

M2_GENES = ["MRC1", "CD163", "MSR1", "ARG1", "TGFB1", "IL10", "VSIG4"]
TREG_GENES = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]
CHECKPT_GENES = ["CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3"]

# Transcript annotation (real_quant products + Ensembl info)
TX_ANNOT = {
    "ENST00000336517.8": "FL 9-exon (canonical)",
    "ENST00000466960.5": "retained-intron",
    "ENST00000394860.3": "5-exon truncated",
    "ENST00000467555.1": "alt-3exon",
    "ENST00000394857.7": "alt-6exon",
    "ENST00000416245.5": "alt-2exon",
    "ENST00000479793.5": "alt UTR",
}


def load_expr(cancer):
    """GBM/LGG 21 gene expression: returns DataFrame index=sample columns=gene (log2TPM)."""
    p = os.path.join(H2, f"TCGA.{cancer}.sampleMap", "HiSeq_TCGA_gene.xena.gz")
    df = pd.read_csv(p, sep="\t", index_col=0, compression="gzip").T
    df.index.name = "sample"
    return df


def score_zs(genes, df):
    """z-score consensus: each gene standardized across samples then averaged. df rows are samples, columns are genes."""
    avail = [g for g in genes if g in df.columns]
    if not avail:
        return pd.Series(np.nan, index=df.index)
    sub = df[avail].astype(float).T       # genes x samples
    v = sub.std(axis=1) > 0
    if not v.any():
        return pd.Series(np.nan, index=df.index)
    z = ((sub.loc[v] - sub.loc[v].mean(axis=1).values[:, None])
         / sub.loc[v].std(axis=1).values[:, None])
    return z.mean(axis=0)


def spearman_p(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan, np.nan, int(m.sum())
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), int(m.sum())


def main():
    print("=== Article 3: ZP3 transcript PSI quantification + immune association ===\n")

    # 1. Load PSI proportion matrix (19131 x 7)
    prop = pd.read_csv(PROP, index_col=0)
    prop = prop[list(TX_ANNOT.keys())]          # Maintain the order of 7 transcripts
    print(f"1. PSI matrix: {prop.shape[0]} samples × {prop.shape[1]} transcripts")

    # Distinguish TCGA tumors (01 segment) and GTEx
    tcga_tumor = [s for s in prop.index
                  if s.startswith("TCGA-") and s.split("-")[3].startswith("01")]
    print(f"   TCGA tumor samples: {len(tcga_tumor)}")

    # 2. Load GBM/LGG expression
    print("\n2. Load GBM/LGG immune expression...")
    expr = {c: load_expr(c) for c in ("GBM", "LGG")}
    for c, df in expr.items():
        print(f"   {c}: {df.shape[0]} samples × {df.shape[1]} genes")

    # == 3. PSI of each transcript in GBM/LGG + comparison between types ==
    print("\n3. GBM vs LGG PSI comparison...")
    psi_status = []
    for c in ("GBM", "LGG"):
        sams = expr[c].index.intersection(prop.index)
        for tx in TX_ANNOT:
            vals = prop.loc[sams, tx].dropna()
            psi_status.append({
                "Cancer": c, "Transcript": tx, "Annotation": TX_ANNOT[tx],
                "N": len(vals), "PSI_median": vals.median(),
                "PSI_mean": vals.mean(), "PSI_p25": vals.quantile(0.25),
                "PSI_p75": vals.quantile(0.75)})
    # MWU GBM vs LGG
    gbm_sams = expr["GBM"].index.intersection(prop.index)
    lgg_sams = expr["LGG"].index.intersection(prop.index)
    for tx in TX_ANNOT:
        g = prop.loc[gbm_sams, tx].dropna().values
        l = prop.loc[lgg_sams, tx].dropna().values
        if len(g) >= 5 and len(l) >= 5:
            u, p = stats.mannwhitneyu(g, l, alternative="two-sided")
            psi_status.append({
                "Cancer": "GBM_vs_LGG", "Transcript": tx,
                "Annotation": TX_ANNOT[tx],
                "N": min(len(g), len(l)),
                "PSI_median": float(np.median(g) - np.median(l)),
                "PSI_mean": float(np.mean(g) - np.mean(l)),
                "MWU_p": p})
    psi_df = pd.DataFrame(psi_status)
    psi_df.to_csv(os.path.join(OUT, "psi_status_by_transcript.csv"), index=False)
    print(f"   Saved psi_status_by_transcript.csv ({len(psi_df)} rows)")

    # == 4. Transcript PSI × immune feature association ==
    print("\n4. Transcript PSI × immune feature association...")
    rec = []
    for c in ("GBM", "LGG"):
        sams = expr[c].index.intersection(prop.index)
        ex = expr[c].loc[sams]
        for gset_name, genes in [("M2", M2_GENES), ("Treg", TREG_GENES),
                                 ("Checkpoint", CHECKPT_GENES)]:
            sc = score_zs(genes, ex)
            for tx in TX_ANNOT:
                x = prop.loc[sams, tx].values
                rho, p, n = spearman_p(x, sc.values)
                rec.append({"Cancer": c, "GeneSet": gset_name,
                            "Transcript": tx, "Annotation": TX_ANNOT[tx],
                            "Rho": rho, "P": p, "N": n})
    corr = pd.DataFrame(rec)
    # BH FDR
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(corr["P"].fillna(1), method="fdr_bh")
    corr["FDR"] = fdr
    corr["Significant"] = corr["FDR"] < 0.05
    corr.to_csv(os.path.join(OUT, "psi_immune_correlation.csv"), index=False)
    print(f"   Saved psi_immune_correlation.csv ({len(corr)} associations)")
    print(f"   Number of significant associations (FDR<0.05): {corr['Significant'].sum()}")
    sig = corr[corr["Significant"]]
    if len(sig):
        print("\n   ------- Significant associations -------")
        print(sig[["Cancer", "GeneSet", "Transcript", "Annotation",
                   "Rho", "P", "FDR", "N"]].to_string(index=False))

    # Key isoform tracking: 5-exon truncated (ENST00000394860.3)
    print("\n=== Key isoform tracking: ENST00000394860.3 (5-exon truncated) ===")
    for c in ("GBM", "LGG"):
        tx = "ENST00000394860.3"
        sams = expr[c].index.intersection(prop.index)
        for gset_name, genes in [("M2", M2_GENES), ("Treg", TREG_GENES),
                                 ("Checkpoint", CHECKPT_GENES)]:
            sc = score_zs(genes, expr[c].loc[sams])
            rho, p, n = spearman_p(prop.loc[sams, tx].values, sc.values)
            print(f"   {c} {gset_name}: ρ={rho:+.3f}, p={p:.3e}, n={n}")

    # == 5. Plotting ==
    print("\n5. Plotting...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # a) GBM vs LGG PSI median comparison (only main 4 transcripts)
    ax = axes[0]
    main_tx = ["ENST00000336517.8", "ENST00000466960.5",
               "ENST00000394860.3", "ENST00000394857.7"]
    gbm_med = [psi_df[(psi_df.Cancer == "GBM") & (psi_df.Transcript == t)]["PSI_median"].values[0]
               if len(psi_df[(psi_df.Cancer == "GBM") & (psi_df.Transcript == t)]) else np.nan
               for t in main_tx]
    lgg_med = [psi_df[(psi_df.Cancer == "LGG") & (psi_df.Transcript == t)]["PSI_median"].values[0]
               if len(psi_df[(psi_df.Cancer == "LGG") & (psi_df.Transcript == t)]) else np.nan
               for t in main_tx]
    x = np.arange(len(main_tx))
    w = 0.35
    ax.bar(x - w/2, gbm_med, w, label="GBM (n=153)", color="#A32D2D", alpha=0.85)
    ax.bar(x + w/2, lgg_med, w, label="LGG (n=509)", color="#1D9E75", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([t.split(".")[0] for t in main_tx], rotation=30, fontsize=8)
    ax.set_ylabel("Median PSI")
    ax.set_title("ZP3 isoform PSI: GBM vs LGG")
    ax.legend()

    # b) PSI × immune correlation heatmap
    ax = axes[1]
    piv = corr.pivot_table(index="Transcript", columns=["Cancer", "GeneSet"],
                           values="Rho")
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([f"{c[0]}-{c[1]}" for c in piv.columns],
                       rotation=45, fontsize=8)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([t.split(".")[0] for t in piv.index], fontsize=8)
    ax.set_title("PSI × immune signature (Spearman ρ)")
    fig.colorbar(im, ax=ax, shrink=0.7)

    # c) Key isoform PSI × immune (mainly 94860.3)
    ax = axes[2]
    tx_focus = "ENST00000394860.3"
    for c, col in [("GBM", "#A32D2D"), ("LGG", "#1D9E75")]:
        sams = expr[c].index.intersection(prop.index)
        psiv = prop.loc[sams, tx_focus].values
        for gset_name in ["M2", "Treg", "Checkpoint"]:
            genes = {"M2": M2_GENES, "Treg": TREG_GENES, "Checkpoint": CHECKPT_GENES}[gset_name]
            sc = score_zs(genes, expr[c].loc[sams])
            ax.scatter(psiv, sc.values, s=6, alpha=0.4, color=col,
                       label=f"{c} {gset_name}" if ax.collections.__len__() < 6 else "")
    ax.set_xlabel(f"PSI {tx_focus.split('.')[0]}")
    ax.set_ylabel("Immune signature (z-score consensus)")
    ax.set_title(f"PSI({tx_focus.split('.')[0]}) vs immune")
    ax.legend(fontsize=7)

    plt.tight_layout()
    fig_path = os.path.join(OUT, "fig_zp3_psi_immune.png")
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    print(f"   Saved {fig_path}")

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
