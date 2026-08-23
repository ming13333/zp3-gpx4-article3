#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ZP3 真实异构体定量（基于 TCGA TARGET GTEx isoform TPM）
=========================================================
数据：本地 TcgaTargetGtex_rsem_isoform_tpm.gz（4.3 GB）
方法：
  1. 流式读取，仅保留 ZP3 候选转录本（Ensembl 正链 8 个转录本）。
  2. 转换 log2(TPM) -> TPM（相对量），计算每样本各转录本占 ZP3 总表达比例。
  3. 用 GDC 映射（tcga_disease_map.json）做 TCGA 癌种分组；
     用 GTEX_phenotype.gz 区分 GTEx 正常样本。
  4. 比较肿瘤 vs 正常中的异构体比例差异（Mann-Whitney U + BH FDR）。

输出：
  - zp3_isoform_proportions.csv：每样本各转录本比例
  - zp3_isoform_tumor_vs_normal.csv：肿瘤 vs 正常差异检验
  - zp3_isoform_by_cancer.csv：各癌种异构体比例中位数
  - fig_zp3_isoform_proportions.png：可视化
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

# ZP3 正链候选转录本（Ensembl GRCh38，从 overlap 查询获得）
ZP3_TRANSCRIPTS = {
    'ENST00000336517', 'ENST00001135277', 'ENST00000394857',
    'ENST00000416245', 'ENST00000394860', 'ENST00000466960',
    'ENST00000479793', 'ENST00000467555'
}

# ---------------------------------------------------------------------------
# 1. 流式读取 isoform TPM（仅 ZP3 候选转录本）
# ---------------------------------------------------------------------------
def read_zp3_isoforms(path, target_prefixes, chunk=5000):
    """流式读取 gz isoform 矩阵，仅保留 target_prefixes（不带版本号）。
    返回 DataFrame：index=转录本 ID（含版本号），columns=样本 id。"""
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
                print(f"    已扫描 {n} 个转录本行...")
    df = pd.DataFrame.from_dict(rows, orient="index", columns=samples)
    return df


# ---------------------------------------------------------------------------
# 2. 加载 GDC 癌种映射 + GTEx 表型
# ---------------------------------------------------------------------------
def load_gdc_map(path):
    with open(path) as f:
        return json.load(f)


def load_gtex_pheno(path):
    """GTEX_phenotype.gz 是 tsv，含 Sample 列和 _primary_site 等。"""
    df = pd.read_csv(path, sep="\t", compression="gzip")
    # 样本列名可能在 'Sample' 或第一列
    sample_col = 'Sample' if 'Sample' in df.columns else df.columns[0]
    df = df.rename(columns={sample_col: 'sample_id'})
    return df


# ---------------------------------------------------------------------------
# 3. 主流程
# ---------------------------------------------------------------------------
def main():
    print("=== ZP3 真实异构体定量（基于 TCGA TARGET GTEx isoform TPM）===\n")
    if not os.path.exists(DATA_ISO):
        print(f"!! 找不到 isoform 数据: {DATA_ISO}"); sys.exit(1)

    # 3.1 读取 isoform
    print("1. 流式读取 isoform TPM（仅 ZP3 候选转录本）...")
    t0 = time.time()
    iso = read_zp3_isoforms(DATA_ISO, ZP3_TRANSCRIPTS)
    print(f"   找到 {iso.shape[0]} 个 ZP3 转录本 × {iso.shape[1]} 个样本，"
          f"耗时 {time.time()-t0:.1f}s")
    if iso.empty:
        print("!! 未找到任何 ZP3 转录本，退出"); sys.exit(1)

    # 3.2 转换 log2(TPM) -> TPM（相对量）
    print("\n2. 转换 log2(TPM) -> TPM 并计算每样本异构体比例...")
    tpm = 2 ** iso  # log2(TPM) 的逆变换，得到相对 TPM
    # 避免负值（理论上 2^log2(TPM) 应为正，但底值 -9.9658 会得到极小值）
    tpm = tpm.clip(lower=0)
    # 每样本总 ZP3 表达
    total_zp3 = tpm.sum(axis=0)
    # 每样本各转录本比例（总 ZP3=0 的样本比例设为 NaN）。
    # 注意：不用 prop.where(cond) —— pandas 会把 Series 条件与行索引对齐
    # 而非列索引，导致全 NaN；改为显式按列赋值。
    prop = tpm.div(total_zp3, axis=1)
    prop.loc[:, total_zp3 <= 0] = np.nan
    print(f"   总 ZP3 TPM=0 的样本数: {(total_zp3 == 0).sum()}")

    # 保存比例矩阵
    prop_path = os.path.join(OUT_DIR, "zp3_isoform_proportions.csv")
    prop.T.to_csv(prop_path)
    print(f"   比例矩阵已保存: {prop_path}")

    # 3.3 样本分类
    samples = list(prop.columns)
    # TCGA 肿瘤：TCGA- 开头且样本段以 01 开头
    tcga_tumor_mask = [s.startswith("TCGA-") and s.split("-")[3].startswith("01")
                       for s in samples]
    tcga_tumor_samples = [s for s, m in zip(samples, tcga_tumor_mask) if m]
    # GTEx 正常
    gtex_samples = [s for s in samples if s.startswith("GTEX-")]
    print(f"\n3. 样本分类:")
    print(f"   TCGA 肿瘤样本: {len(tcga_tumor_samples)}")
    print(f"   GTEx 正常样本: {len(gtex_samples)}")

    # 3.4 肿瘤 vs 正常：各转录本比例差异
    print("\n4. 肿瘤 vs 正常（GTEx）异构体比例差异...")
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
    print(f"   结果已保存: {res_path}")
    print(res_df.to_string(index=False))

    # 3.5 按癌种分析（TCGA 肿瘤）
    print("\n5. 按癌种分析异构体比例...")
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
        print(f"   覆盖 {len(cancers)} 个癌种")

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
        print(f"   结果已保存: {cancer_path}")
    else:
        print(f"   (跳过：GDC 映射缓存不存在: {GDC_MAP})")

    # 3.6 可视化
    print("\n6. 生成可视化...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # (a) 肿瘤 vs 正常箱线图
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

    # (b) 转录本表达量（log2 TPM）小提琴图
    log_data = []
    for tid in iso.index:
        for s in tcga_tumor_samples[:500]:  # 采样避免图过密
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
    print(f"   图表已保存: {fig_path}")

    print("\n=== 分析完成 ===")


if __name__ == "__main__":
    main()
