#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Article 3：ZP3 转录本 PSI 定量 + 免疫关联
============================================
基于已就位的 Toil isoform TPM（TcgaTargetGtex_rsem_isoform_tpm.gz, 4.2G）
对 ZP3 7 个转录本做 PSI（Percent Spliced In）定量：
  PSI_t = TPM(transcript t) / Σ_all ZP3 transcripts TPM

分析：
  1. 全 TCGA 肿瘤 vs GTEx 正常：PSI 差异（已有 tumor_vs_normal.csv，此处聚焦 GBM/LGG）
  2. GBM(153) vs LGG(509)：PSI 类型间比较
  3. **转录本级别免疫关联**：GBM/LGG 中每转录本 PSI 与免疫特征
     (M2 / Treg / Checkpoint z-score 共识) 的 Spearman 相关 —— 回答
     "哪个转录本驱动 ZP3-免疫抑制关联"
  4. 关键异构体：ENST00000394860.3（5 外显子截短、肿瘤富集 22 倍）的
     PSI 特异性关联

数据：
  - zp3_isoform_proportions.csv（19131 样本 × 7 转录本，已由 real_quant 产出）
  - h2_bulk/TCGA.GBM.sampleMap 与 TCGA.LGG.sampleMap（21 免疫基因 log2TPM）

产物：
  psi_status_by_transcript.csv  —— GBM/LGG 各转录本 PSI 中位数 + GBM vs LGG MWU
  psi_immune_correlation.csv    —— 各转录本 PSI × 免疫特征 Spearman
  fig_zp3_psi_immune.png        —— 可视化
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

# 转录本注释（real_quant 产物 + Ensembl 信息）
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
    """GBM/LGG 21 基因表达：返回 DataFrame index=sample columns=gene (log2TPM)。"""
    p = os.path.join(H2, f"TCGA.{cancer}.sampleMap", "HiSeq_TCGA_gene.xena.gz")
    df = pd.read_csv(p, sep="\t", index_col=0, compression="gzip").T
    df.index.name = "sample"
    return df


def score_zs(genes, df):
    """z-score 共识：每基因跨样本标准化后取均值。df 行为样本、列为基因。"""
    avail = [g for g in genes if g in df.columns]
    if not avail:
        return pd.Series(np.nan, index=df.index)
    sub = df[avail].astype(float).T       # 基因×样本
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
    print("=== Article 3: ZP3 转录本 PSI 定量 + 免疫关联 ===\n")

    # 1. 载入 PSI 比例矩阵（19131 × 7）
    prop = pd.read_csv(PROP, index_col=0)
    prop = prop[list(TX_ANNOT.keys())]          # 保持 7 转录本顺序
    print(f"1. PSI 矩阵: {prop.shape[0]} 样本 × {prop.shape[1]} 转录本")

    # 判别 TCGA 肿瘤（01 段）与 GTEx
    tcga_tumor = [s for s in prop.index
                  if s.startswith("TCGA-") and s.split("-")[3].startswith("01")]
    print(f"   TCGA 肿瘤样本: {len(tcga_tumor)}")

    # 2. 加载 GBM/LGG 表达
    print("\n2. 加载 GBM/LGG 免疫表达...")
    expr = {c: load_expr(c) for c in ("GBM", "LGG")}
    for c, df in expr.items():
        print(f"   {c}: {df.shape[0]} 样本 × {df.shape[1]} 基因")

    # == 3. GBM/LGG 各转录本 PSI + 类型间比较 ==
    print("\n3. GBM vs LGG PSI 比较...")
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
    print(f"   已保存 psi_status_by_transcript.csv ({len(psi_df)} 行)")

    # == 4. 转录本 PSI × 免疫特征关联 ==
    print("\n4. 转录本 PSI × 免疫特征关联...")
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
    print(f"   已保存 psi_immune_correlation.csv ({len(corr)} 条关联)")
    print(f"   显著关联数(FDR<0.05): {corr['Significant'].sum()}")
    sig = corr[corr["Significant"]]
    if len(sig):
        print("\n   ------- 显著关联 -------")
        print(sig[["Cancer", "GeneSet", "Transcript", "Annotation",
                   "Rho", "P", "FDR", "N"]].to_string(index=False))

    # 关键异构体追踪：5-exon truncated (ENST00000394860.3)
    print("\n=== 关键异构体追踪: ENST00000394860.3 (5-exon truncated) ===")
    for c in ("GBM", "LGG"):
        tx = "ENST00000394860.3"
        sams = expr[c].index.intersection(prop.index)
        for gset_name, genes in [("M2", M2_GENES), ("Treg", TREG_GENES),
                                 ("Checkpoint", CHECKPT_GENES)]:
            sc = score_zs(genes, expr[c].loc[sams])
            rho, p, n = spearman_p(prop.loc[sams, tx].values, sc.values)
            print(f"   {c} {gset_name}: ρ={rho:+.3f}, p={p:.3e}, n={n}")

    # == 5. 绘图 ==
    print("\n5. 绘图...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # a) GBM vs LGG PSI 中位对比（仅主 4 转录本）
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

    # b) PSI×免疫关联热图
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

    # c) 关键异构体 PSI×免疫（主要看 94860.3）
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
    print(f"   已保存 {fig_path}")

    print("\n=== 完成 ===")


if __name__ == "__main__":
    main()