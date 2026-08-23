# -*- coding: utf-8 -*-
"""
Article 3 — 跨癌种 ZP3 异构体 PSI 特异性比较
=============================================
在 TCGA 泛癌中比较 7 个 ZP3 转录本的 PSI（percent spliced in）指纹，并检验
转录本 PSI 与免疫特征关联的跨癌种模式：

  1) PSI 指纹热图（32 癌种 × 7 转录本，中位 PSI）
  2) per-cancer 转录本 PSI × 免疫特征 Spearman（32 × 7 × 7 = 1568 条）
  3) 生态学分析：癌种级 FL-PSI / RI-PSI 中位 × ZP3-免疫关联强度(Avg_Rho)
  4) 异构体切换指数 log2(FL_PSI / RI_PSI) 跨癌种排序

数据源（全部本地）：
  - zp3_isoform_proportions.csv （19131 样本 × 7 转录本 PSI）
  - TcgaTargetGtex_rsem_gene_tpm.gz（1.3G 本地 TPM，免疫基因提取）
  - tcga_disease_map.json（样本 barcode -> 癌种）
  - ensg_map.json（symbol -> Ensembl 缓存）

产物目录：zp3_psi_pancancer_results/
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
ROOT = os.path.dirname(os.path.dirname(BASE))      # 项目根
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_pancancer_results")
os.makedirs(OUT, exist_ok=True)

DATA_TPM = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "phase1_knowledge_gap_filling", "TcgaTargetGtex_rsem_gene_tpm.gz")
PROP_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_isoform_proportions.csv")
DISEASE_MAP = os.path.join(ROOT, "output", "tcga_pancan", "tcga_disease_map.json")
ENSG_CACHE = os.path.join(ROOT, "output", "tcga_pancan", "ensg_map.json")
SUMMARY_CSV = os.path.join(ROOT, "output", "tcga_pancan", "tcga_pancan_cancer_summary.csv")

# 与 tcga_pancan_zp3_analysis.py 完全一致的免疫基因集
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
# 转录本语义标签
TX_LABEL = {
    "ENST00000336517.8": "FL canonical (9-exon)",
    "ENST00000466960.5": "Retained-intron",
    "ENST00000394860.3": "5-exon truncated",
    "ENST00000394857.7": "Alt-exon isoform",
    "ENST00000467555.1": "Short isoform",
    "ENST00000416245.5": "Mid isoform",
    "ENST00000479793.5": "Alt-terminal isoform",
}
FL = "ENST00000336517.8"   # 经典全长
RI = "ENST00000466960.5"   # retained-intron


def get_ensg_map(symbols):
    with open(ENSG_CACHE) as f:
        cache = json.load(f)
    return {s: cache.get(s) for s in symbols}


def read_target_genes(path, target_strip):
    """流式读取 TPM 中目标基因（target_strip: ensg去版本 -> symbol 或 ensg）。"""
    rows = {}
    with gzip.open(path, "rt") as f:
        first = f.readline()
        samples = first.rstrip("\n").split("\t")[1:]
        print(f"    共 {len(samples)} 个样本")
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
                print(f"    已扫描 {n} 个基因行...")
    return pd.DataFrame.from_dict(rows, orient="index", columns=samples)


def zscore_consensus_score(mat_gx_s, ensgs):
    """z-score 共识：每基因跨样本标准化后取均值（与泛癌脚本一致）。返回样本级 Series。"""
    sub = mat_gx_s.loc[ensgs]
    gene_mean = sub.mean(axis=1)
    gene_std = sub.std(axis=1)
    valid = gene_std > 0
    if not valid.any():
        return None
    z = (sub.loc[valid] - gene_mean[valid].values[:, None]) / gene_std[valid].values[:, None]
    return z.mean(axis=0)


def main():
    print("=== Article 3: 跨癌种 ZP3 异构体 PSI 特异性比较 ===\n")

    # ---- 1. PSI 矩阵 ----
    print("1. 读取 PSI 比例矩阵 ...")
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    print(f"   PSI 矩阵 {psi.shape[0]} 样本 × {psi.shape[1]} 转录本")

    # ---- 2. 免疫基因 TPM（流式）----
    all_symbols = sorted({g for s in IMMUNE_GENE_SETS.values() for g in s})
    sym2ensg = get_ensg_map(all_symbols)
    unresolved = [s for s, e in sym2ensg.items() if not e]
    if unresolved:
        print(f"  !! 未解析免疫基因: {unresolved}")
    target_ids = [e for e in sym2ensg.values() if e]
    print(f"2. 流式读取真实 TPM（{len(target_ids)} 个免疫基因）...")
    t0 = time.time()
    mat = read_target_genes(DATA_TPM, {e: e for e in target_ids})
    print(f"   读取完成 {mat.shape[0]} 基因 × {mat.shape[1]} 样本，耗时 {time.time()-t0:.1f}s")

    # ---- 3. 免疫评分（全 TCGA 肿瘤样本）----
    samples = list(mat.columns)
    tcga_mask = [s.startswith("TCGA-") and s.split("-")[3].startswith("01") for s in samples]
    tcga_samples = [s for s, m in zip(samples, tcga_mask) if m]
    mat_t = mat[tcga_samples]
    print(f"3. TCGA 肿瘤样本 {len(tcga_samples)} 个，计算 7 特征 z-score 共识评分 ...")

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
    print(f"   完成评分，{len(score_mat)} 个样本有评分")

    # ---- 4. join PSI + 免疫评分 + 癌种 ----
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
    print(f"4. 合并后 {len(df)} 样本（PSI × 免疫评分 × 癌种）")
    if len(df) == 0:
        print("!! 合并为空，退出"); sys.exit(1)
    df.to_csv(os.path.join(OUT, "psi_immune_joined_samples.csv"), index=False)

    # ---- 5. per-cancer 转录本 PSI × 免疫特征关联 ----
    print("5. per-cancer 转录本 PSI × 免疫特征 Spearman ...")
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
    print(f"   {len(corr)} 条关联，覆盖 {corr['Cancer'].nunique()} 癌种")
    n_sig = (corr["FDR"] < 0.05).sum()
    print(f"   显著(FDR<0.05): {n_sig} 条")

    # FL 转录本跨癌种方向一致性
    fl = corr[corr["Transcript"] == FL]
    fl_pos_sig = ((fl["FDR"] < 0.05) & (fl["Rho"] > 0)).sum()
    fl_neg_sig = ((fl["FDR"] < 0.05) & (fl["Rho"] < 0)).sum()
    fl_all = len(fl)
    print(f"   FL 关联: {fl_all} 条，正显著 {fl_pos_sig} / 负显著 {fl_neg_sig}")
    ri = corr[corr["Transcript"] == RI]
    ri_pos_sig = ((ri["FDR"] < 0.05) & (ri["Rho"] > 0)).sum()
    ri_neg_sig = ((ri["FDR"] < 0.05) & (ri["Rho"] < 0)).sum()
    print(f"   RI 关联: {len(ri)} 条，正显著 {ri_pos_sig} / 负显著 {ri_neg_sig}")

    # ---- 6. 癌种级 PSI 指纹 + 生态学 ----
    print("6. 癌种级 PSI 指纹与生态学分析 ...")
    fp = df.groupby("Cancer")[list(psi.columns)].median()
    fp["N"] = df.groupby("Cancer").size()
    fp = fp.reset_index().rename(columns={"index": "Cancer"})
    fp = fp.sort_values(FL, ascending=False)
    fp.to_csv(os.path.join(OUT, "psi_pancancer_fingerprint.csv"), index=False)

    # 切换指数 per-sample 后取癌种中位
    df["switch_index"] = np.log2((df[FL].clip(lower=1e-6)) / (df[RI].clip(lower=1e-6)))
    sw = df.groupby("Cancer")["switch_index"].median().reset_index()

    # 生态学：与泛癌 Avg_Rho 合并
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
            print(f"   生态学相关: {lab}-PSI 中位 × ZP3-免疫 Avg_Rho: ρ={rho:+.3f}, p={p:.3f}, n={m.sum()}")
        x2 = eco["switch_index"].values
        rho2, p2 = stats.spearmanr(x2[m], y[m])
        print(f"   生态学相关: switch_index × Avg_Rho: ρ={rho2:+.3f}, p={p2:.3f}")

    # ---- 7. 图 ----
    print("7. 绘图 ...")
    plot_figure(fp, eco, corr)
    print("\n=== 完成，产物在", OUT, "===")


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

    # (a) PSI 指纹热图
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

    # (b) 生态学：FL-PSI × Avg_Rho
    ax = fig.add_subplot(2, 2, 2)
    x = eco[FL].values.astype(float); y = eco["Avg_Rho"].values.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    rho, p = stats.spearmanr(x[m], y[m])
    ax.scatter(x[m], y[m], s=55, alpha=0.8, c="#378ADD", edgecolor="white")
    # 高亮胶质瘤
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

    # (c) 生态学：RI-PSI × Avg_Rho
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

    # (d) FL 转录本 PSI × 免疫特征 per-cancer 关联热图
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
