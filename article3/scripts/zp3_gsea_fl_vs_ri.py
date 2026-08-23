# -*- coding: utf-8 -*-
"""
Article 3 补强② — FL-high vs RI-high 的 GSEA 功能注释
=====================================================
目的：检验异构体切换的功能暗示——FL 全长占比高的样本 vs RI (retained-intron)
占比高的样本，差异表达基因富集哪些通路（重点：免疫抑制通路）。

设计：
  1. 样本：GBM+LGG 合并胶质瘤（TCGA），按 PSI 分组
     - FL_high: FL-PSI ≥ 上三分位 且 RI-PSI ≤ 下三分位（纯 FL 组）
     - RI_high: RI-PSI ≥ 上三分位 且 FL-PSI ≤ 下三分位（纯 RI 组）
  2. 差异表达：两组间逐基因 Wilcoxon 秩和检验 → z 得分作为排序指标
     （z = 正态近似统计量，正值 = FL_high 中上调）
  3. GSEA: gseapy.prerank，基因集 MSigDB Hallmark_2020 + C7 免疫
  4. 输出：富集表 + 关键免疫抑制通路的富集图

产物：zp3_gsea_results/
"""
import os
import json
import gzip
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_gsea_results")
os.makedirs(OUT, exist_ok=True)

DATA_TPM = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "phase1_knowledge_gap_filling", "TcgaTargetGtex_rsem_gene_tpm.gz")
PROP_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_isoform_proportions.csv")
DISEASE_MAP = os.path.join(os.path.dirname(os.path.dirname(BASE)),
                           "output", "tcga_pancan", "tcga_disease_map.json")

FL = "ENST00000336517.8"   # 经典全长
RI = "ENST00000466960.5"   # retained-intron
HALLMARK = "MSigDB_Hallmark_2020"
C7_IMMUNE = "MSigDB_Immunologic_Signatures"   # 若不可用则跳过
ENSG_SYMBOL_CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "phase1_knowledge_gap_filling", "ensg_symbol_map.json")


def read_psi():
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    return psi


def load_expression(target_ensg_set, sample_subset):
    """流式读取 TPM 中目标基因 × 目标样本（返回 genes × samples）。"""
    rows = {}
    sample_list = list(sample_subset)
    sample_set = set(sample_list)
    with gzip.open(DATA_TPM, "rt") as f:
        first = f.readline()
        all_samples = first.rstrip("\n").split("\t")[1:]
        idx = [i for i, s in enumerate(all_samples) if s in sample_set]
        n = 0
        while True:
            lines = f.readlines(65536)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip("\n").split("\t")
                gid = parts[0].split(".")[0]
                if gid in target_ensg_set:
                    rows[gid] = [float(parts[i + 1]) for i in idx]
            n += len(lines)
            if n % 400000 == 0:
                print(f"    已扫描 {n} 行...")
    return pd.DataFrame(rows, index=sample_list).T  # genes × samples


def main():
    print("=== Article 3 补强②: FL-high vs RI-high GSEA ===\n")

    # 1. PSI + 癌种
    psi = read_psi()
    with open(DISEASE_MAP) as f:
        p2c = json.load(f)
    cancer_of = {}
    for s in psi.index:
        if s.startswith("TCGA-") and len(s.split("-")) >= 4:
            cancer_of[s] = p2c.get("-".join(s.split("-")[:3]), "UNKNOWN")
        else:
            cancer_of[s] = "UNKNOWN"

    gl = [s for s, c in cancer_of.items() if c in ("GBM", "LGG")]
    print(f"胶质瘤样本 (GBM+LGG): {len(gl)}")

    sub = psi.loc[gl, [FL, RI]].dropna()
    print(f"有 PSI 值: {len(sub)}")

    # 2. 三分位分组（纯组）
    fl_q = sub[FL].quantile(0.67)
    ri_q = sub[RI].quantile(0.67)
    fl_low_q = sub[FL].quantile(0.33)
    ri_low_q = sub[RI].quantile(0.33)
    print(f"FL 三分位: {sub[FL].quantile(0.33):.3f} / {fl_q:.3f}")
    print(f"RI 三分位: {sub[RI].quantile(0.33):.3f} / {ri_q:.3f}")

    fl_high = sub[(sub[FL] >= fl_q) & (sub[RI] <= ri_low_q)].index
    ri_high = sub[(sub[RI] >= ri_q) & (sub[FL] <= fl_low_q)].index
    print(f"FL_high (纯 FL): {len(fl_high)} | RI_high (纯 RI): {len(ri_high)}")
    if len(fl_high) < 15 or len(ri_high) < 15:
        print("!! 组样本不足，放宽条件")
        fl_high = sub[sub[FL] >= fl_q].index
        ri_high = sub[sub[RI] >= ri_q].index
        print(f"  放宽后: FL_high={len(fl_high)} RI_high={len(ri_high)}")

    # 3. 读全表达（两组样本）
    print("读取 TPM（两组样本）...")
    all_samp = sorted(set(fl_high) | set(ri_high))
    # 目标基因集 = 免疫相关（加快读取；GSEA 需要完整基因空间时再扩）
    # 这里先读全部基因（~60000 行，流式可行）
    ensg_target = None  # None = 读所有基因
    # 但读所有基因内存大，折中：跳过表达过低的基因行会降低质量，全部读
    rows = {}
    idx_map = {}
    with gzip.open(DATA_TPM, "rt") as f:
        first = f.readline()
        all_samples = first.rstrip("\n").split("\t")[1:]
        idx = [i for i, s in enumerate(all_samples) if s in set(all_samp)]
        sample_idx = {all_samples[i]: i for i in idx}
        n = 0
        while True:
            lines = f.readlines(65536)
            if not lines:
                break
            for ln in lines:
                parts = ln.rstrip("\n").split("\t")
                gid = parts[0].split(".")[0]
                rows[gid] = [float(parts[i + 1]) for i in idx]
            n += len(lines)
            if n % 400000 == 0:
                print(f"    已扫描 {n} 行...")
    expr = pd.DataFrame(rows, index=all_samp).T  # genes × samples
    print(f"表达矩阵: {expr.shape[0]} 基因 × {expr.shape[1]} 样本")

    # 4. 差异表达（Wilcoxon，两组）
    print("差异表达检验（Wilcoxon）...")
    g_fl = expr[fl_high]
    g_ri = expr[ri_high]
    genes = list(expr.index)
    recs = []
    for g in genes:
        a = g_fl.loc[g].values.astype(float)
        b = g_ri.loc[g].values.astype(float)
        # TPM 文件是 log2(TPM+ε)，过滤极低表达（log2 均值 < -3 ≈ TPM<0.12）
        if np.nanmean(a) < -3.0 and np.nanmean(b) < -3.0:
            continue  # 过滤低表达
        try:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            # z 近似
            na, nb = len(a), len(b)
            mu = na * nb / 2
            sd = np.sqrt(na * nb * (na + nb + 1) / 12)
            z = (u - mu) / sd if sd > 0 else 0.0
            recs.append({"gene": g, "z": z, "p": p,
                         "mean_FL": float(np.nanmean(a)), "mean_RI": float(np.nanmean(b))})
        except Exception:
            continue
    de = pd.DataFrame(recs)
    print(f"DE 基因数: {len(de)}")

    # 基因名：ensg -> symbol。优先用本地全基因组映射缓存；无则退回 ensg_map。
    ensg2sym = {}
    if os.path.exists(ENSG_SYMBOL_CACHE):
        with open(ENSG_SYMBOL_CACHE) as f:
            ensg2sym = json.load(f)
        print(f"加载全基因组 symbol 映射: {len(ensg2sym)} 条")
    else:
        print("!! 无全基因组映射缓存，退回小映射（仅覆盖免疫基因）")
        ensg_map_path = os.path.join(os.path.dirname(os.path.dirname(BASE)),
                                     "output", "tcga_pancan", "ensg_map.json")
        with open(ensg_map_path) as f:
            sym2ensg = json.load(f)
        ensg2sym = {v: k for k, v in sym2ensg.items() if v}
    de["symbol"] = de["gene"].map(lambda g: ensg2sym.get(g, g))
    n_mapped = (de["symbol"] != de["gene"]).sum()
    print(f"symbol 映射: {n_mapped}/{len(de)} 基因成功 ({n_mapped/len(de)*100:.1f}%)")

    # 5. prerank GSEA
    rnk = de[["symbol", "z"]].drop_duplicates("symbol").sort_values("z", ascending=False)
    rnk = rnk.dropna()
    print(f"排序列表: {len(rnk)} 基因")
    rnk.to_csv(os.path.join(OUT, "fl_vs_ri_ranklist.csv"), index=False)

    gene_sets = [HALLMARK]
    try:
        gp.get_library(C7_IMMUNE, organism="Human")
        gene_sets.append(C7_IMMUNE)
    except Exception:
        print(f"  {C7_IMMUNE} 不可用，仅用 Hallmark")

    print("运行 prerank GSEA ...")
    res = gp.prerank(rnk=rnk, gene_sets=gene_sets,
                     outdir=os.path.join(OUT, "gsea_out"),
                     min_size=5, max_size=500, permutation_num=1000,
                     seed=42, threads=4, no_plot=True)
    # 汇总结果（gseapy 1.3.1 输出文件名：gseapy.gene_set.prerank.report.csv）
    report = os.path.join(OUT, "gsea_out", "gseapy.gene_set.prerank.report.csv")
    if os.path.exists(report):
        t = pd.read_csv(report)
        t.to_csv(os.path.join(OUT, "gsea_summary.csv"), index=False)
        # Term 列含 "库名__通路名"，拆出通路名
        t["Pathway"] = t["Term"].str.split("__", expand=True)[1].fillna(t["Term"])
        sig = t[t["FWER p-val"] < 0.25].sort_values("NES", ascending=False)
        print(f"\n富集结果: {len(t)} 条 | FWER<0.05: {(t['FWER p-val'] < 0.05).sum()} 条 | FWER<0.25: {len(sig)} 条")
        print("\n=== FL_high 上调通路 (NES>0) ===")
        up = sig[sig["NES"] > 0][["Pathway", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]]
        print(up.head(12).to_string(index=False))
        print("\n=== FL_high 下调通路 (NES<0, RI_high 富集) ===")
        dn = sig[sig["NES"] < 0][["Pathway", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]]
        print(dn.head(12).to_string(index=False))
    else:
        print(f"!! 未找到 GSEA 报告: {report}")
        print("   实际输出:", os.listdir(os.path.join(OUT, "gsea_out")) if os.path.exists(os.path.join(OUT, "gsea_out")) else "无目录")

    print(f"\n产物目录: {OUT}")

    # ---- 图：显著通路 NES 条形图 ----
    if os.path.exists(os.path.join(OUT, "gsea_summary.csv")):
        t = pd.read_csv(os.path.join(OUT, "gsea_summary.csv"))
        t["Pathway"] = t["Term"].str.split("__", expand=True)[1].fillna(t["Term"])
        sig = t[t["FWER p-val"] < 0.25].copy()
        if len(sig):
            sig = sig.sort_values("NES")
            fig, ax = plt.subplots(figsize=(9, max(6, len(sig) * 0.35)))
            colors = ["#C00000" if v > 0 else "#1D9E75" for v in sig["NES"]]
            ax.barh(sig["Pathway"], sig["NES"], color=colors, alpha=0.85,
                    edgecolor="white", linewidth=0.5)
            ax.axvline(0, color="grey", lw=0.8, ls="--")
            ax.set_xlabel("Normalized Enrichment Score (FL_high vs RI_high)", fontsize=11)
            ax.set_title("GSEA Hallmark: FL-high vs RI-high (glioma, n=366)\n"
                         "red = enriched in FL_high | green = enriched in RI_high",
                         fontsize=12)
            ax.grid(axis="x", alpha=0.3)
            fig.tight_layout()
            fig.savefig(os.path.join(OUT, "fig_gsea_fl_vs_ri.png"), dpi=200,
                        bbox_inches="tight")
            plt.close(fig)
            print(f"图已存: fig_gsea_fl_vs_ri.png ({len(sig)} 条显著通路)")


if __name__ == "__main__":
    main()
