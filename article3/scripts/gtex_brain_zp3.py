#!/usr/bin/env python3
"""
GTEx 正常脑对照分析: ZP3 正常 vs 肿瘤表达差异

数据源:
1. GTEX_phenotype (Xena Toil hub) — GTEx 样本组织映射
2. Toil isoform TPM (下载中) — ZP3 基因总 TPM = sum(isoform TPM)
3. cBioPortal GBM/LGG gene expression (已有)
4. HPA protein expression (已有)

产出: ZP3 正常脑 vs GBM vs LGG 表达对比 + 差异倍数
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys, time, gzip, io, requests, json, zlib

DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "phase1_knowledge_gap_filling")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article3", "results", "gtex_normal_brain")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================
# 1. 下载 GTEX phenotype
# ======================
def download_gtex_phenotype():
    """下载 GTEx 表型数据 (样本→组织)"""
    print("[1/5] 下载 GTEx 表型数据...")
    url = "https://toil.xenahubs.net/download/GTEX_phenotype.gz"
    local_path = os.path.join(DIR, "GTEX_phenotype.tsv")
    
    if os.path.exists(local_path) and os.path.getsize(local_path) > 1000:
        print(f"  已有: {local_path}")
        return local_path
    
    print("  流式下载 + 解压...")
    r = requests.get(url, stream=True, timeout=(30, 120))
    r.raise_for_status()
    
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buffer = b""
    text_lines = []
    
    for chunk in r.iter_content(chunk_size=65536):
        buffer += d.decompress(chunk)
        while b"\n" in buffer:
            line, buffer = buffer.split(b"\n", 1)
            text_lines.append(line.decode("utf-8", errors="replace"))
    
    # 写入文件
    with open(local_path, "w", encoding="utf-8") as f:
        for line in text_lines:
            f.write(line + "\n")
    
    print(f"  已保存: {local_path} ({len(text_lines)} 行)")
    return local_path

def load_phenotype(path):
    """加载 GTEx 表型, 筛选脑组织样本"""
    print("\n[2/5] 筛选正常脑组织样本...")
    df = pd.read_csv(path, sep="\t", index_col=0)
    print(f"  总 GTEx 样本: {len(df)}")
    print(f"  列: {list(df.columns)}")
    
    # 找组织相关列
    tissue_cols = [c for c in df.columns if "tissue" in c.lower() or "site" in c.lower() or "body" in c.lower()]
    if tissue_cols:
        print(f"  组织相关列: {tissue_cols[:5]}")
    
    # 找脑组织
    brain_col = None
    for col in df.columns:
        # 查看该列是否含 brain
        if df[col].dtype == object:
            brain_hits = df[col].astype(str).str.lower().str.contains("brain|brain_|cerebr|hippocamp|cortex|frontal|temporal|cerebell", na=False)
            if brain_hits.sum() > 10:
                brain_col = col
                break
    
    if brain_col:
        brain_terms = df[brain_col].astype(str).str.lower()
        is_brain = brain_terms.str.contains("brain|brain_|cerebr|hippocamp|cortex|frontal|temporal|cerebell", na=False)
        brain_samples = df[is_brain]
        print(f"\n  脑组织列: {brain_col}")
        print(f"  脑组织样本: {len(brain_samples)}")
        
        # 细分脑区
        if df[brain_col].dtype == object:
            brain_regions = brain_samples[brain_col].value_counts().head(15)
            print(f"  脑区分部:\n{brain_regions.to_string()}")
        
        return brain_samples, brain_col
    else:
        # 尝试所有列
        print("  未找到专用组织列, 打印前几行样本...")
        print(f"  前5行: {df.head(3)}")
        return df, None

# ======================
# 2. 从已下载的 ZP3 isoform 数据提取 GTEx 样本表达
# ======================
def extract_gtex_zp3_expression(isoform_tsv, brain_samples):
    """从 Toil isoform TPM 提取 GTEx 脑样本的 ZP3 总表达"""
    print("\n[3/5] 提取 GTEx 脑样本 ZP3 表达...")
    
    isoform_path = os.path.join(DIR, isoform_tsv)
    if not os.path.exists(isoform_path):
        print(f"  ⚠ isoform 数据尚未下载: {isoform_path}")
        return None
    
    df_iso = pd.read_csv(isoform_path, sep="\t", index_col=0)
    df_iso = df_iso.T  # 样本 × 转录本
    
    # 筛选 GTEx 样本
    gtex_idx = [s for s in df_iso.index if str(s).startswith("GTEX")]
    df_gtex = df_iso.loc[gtex_idx]
    
    # 计算总 TPM (log2 scale → TPM)
    for col in df_gtex.columns:
        df_gtex[f"{col}_tpm"] = 2**df_gtex[col] - 0.001
    
    tpm_cols = [c for c in df_gtex.columns if c.endswith("_tpm")]
    df_gtex["ZP3_TPM"] = df_gtex[tpm_cols].sum(axis=1)
    df_gtex["ZP3_TPM_log2"] = np.log2(df_gtex["ZP3_TPM"] + 0.001)
    
    # 匹配脑组织样本
    gtex_ids = set(brain_samples.index) if brain_samples is not None else set()
    brain_in_gtex = df_gtex.index.intersection(gtex_ids)
    print(f"  脑组织样本在 isoform 数据中: {len(brain_in_gtex)}/{len(gtex_ids)}")
    
    df_brain = df_gtex.loc[brain_in_gtex] if len(brain_in_gtex) > 0 else df_gtex
    
    return df_brain, df_gtex

# ======================
# 3. 加载 cBioPortal GBM/LGG ZP3 表达
# ======================
def load_tcga_zp3_expression():
    """从已有 cBioPortal 数据加载 TCGA GBM/LGG ZP3"""
    print("\n[4/5] 加载 TCGA GBM/LGG ZP3 表达 (cBioPortal)...")
    
    # 从已有文件读取或在新建
    # 使用 cBioPortal API 获取
    results = {}
    
    for study, entrez_id in [("gbm_tcga", 8277), ("lgg_tcga", 8277)]:
        try:
            import requests as rq
            url = f"https://www.cbioportal.org/api/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data"
            params = {
                'entrezGeneId': entrez_id,
                'sampleListId': f'{study}_rna_seq_v2_mrna'
            }
            r = rq.get(url, params=params, timeout=60, headers={'Accept': 'application/json'})
            if r.status_code == 200:
                data = r.json()
                values = {d['sampleId']: float(d['value']) for d in data}
                results[study] = values
                print(f"  {study}: {len(values)} 样本")
        except Exception as e:
            print(f"  {study}: 失败 {e}")
    
    # 合并
    df_tcga = pd.DataFrame()
    for study, values in results.items():
        proj = "GBM" if "gbm" in study else "LGG"
        s = pd.Series(values, name=f"ZP3_{proj}")
        df_tcga = pd.concat([df_tcga, s.to_frame()], axis=1)
    
    if "ZP3_GBM" in df_tcga.columns and "ZP3_LGG" in df_tcga.columns:
        df_tcga["ZP3_TPM"] = df_tcga["ZP3_GBM"].fillna(df_tcga["ZP3_LGG"])
        df_tcga["Project"] = df_tcga.apply(lambda r: "GBM" if not pd.isna(r.get("ZP3_GBM")) else "LGG", axis=1)
    elif "ZP3_GBM" in df_tcga.columns:
        df_tcga["ZP3_TPM"] = df_tcga["ZP3_GBM"]
        df_tcga["Project"] = "GBM"
    elif "ZP3_LGG" in df_tcga.columns:
        df_tcga["ZP3_TPM"] = df_tcga["ZP3_LGG"]
        df_tcga["Project"] = "LGG"
    
    # log2
    df_tcga["ZP3_TPM_log2"] = np.log2(df_tcga["ZP3_TPM"] + 0.001)
    
    print(f"  TCGA 总计: {len(df_tcga)} 样本")
    return df_tcga

# ======================
# 4. 整合 + 可视化
# ======================
def create_comparison(df_brain, df_tcga, gtex_all):
    """创建正常 vs 肿瘤对比可视化 + 统计"""
    print("\n[5/5] 正常 vs 肿瘤对比分析...")
    
    fig = plt.figure(figsize=(14, 12))
    gs = plt.matplotlib.gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # A: 箱线图 — GTEx Brain vs GBM vs LGG
    ax1 = fig.add_subplot(gs[0, :])
    
    plot_data = []
    labels = []
    colors = []
    
    # GTEx 脑
    if df_brain is not None and len(df_brain) > 0:
        vals = df_brain["ZP3_TPM" if "ZP3_TPM" in df_brain.columns else df_brain.columns[0]]
        plot_data.append(vals.dropna().values)
        labels.append(f"GTEx Brain\n(n={len(vals.dropna())})")
        colors.append("#2ECC71")
    
    # GBM
    if df_tcga is not None and "ZP3_TPM" in df_tcga.columns:
        gbm_vals = df_tcga[df_tcga["Project"] == "GBM"]["ZP3_TPM"].dropna()
        if len(gbm_vals) > 0:
            plot_data.append(gbm_vals.values)
            labels.append(f"TCGA-GBM\n(n={len(gbm_vals)})")
            colors.append("#E74C3C")
    
    # LGG
    if df_tcga is not None and "ZP3_TPM" in df_tcga.columns:
        lgg_vals = df_tcga[df_tcga["Project"] == "LGG"]["ZP3_TPM"].dropna()
        if len(lgg_vals) > 0:
            plot_data.append(lgg_vals.values)
            labels.append(f"TCGA-LGG\n(n={len(lgg_vals)})")
            colors.append("#3498DB")
    
    # 箱线图
    bp = ax1.boxplot(plot_data, patch_artist=True, widths=0.5, showfliers=True, flierprops={'alpha': 0.3, 'markersize': 3})
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    
    ax1.set_xticklabels(labels, fontsize=12)
    ax1.set_ylabel("ZP3 Expression (TPM/FPM)", fontsize=13)
    ax1.set_title("A. ZP3 Expression: Normal Brain vs GBM vs LGG", fontsize=14, fontweight='bold')
    ax1.set_yscale("log")
    
    # 计算差异倍数
    if len(plot_data) >= 2:  # GTEx + GBM
        gtex_median = np.median(plot_data[0]) if len(plot_data) > 0 else None
        fold_changes = []
        for i in range(1, len(plot_data)):
            fc = np.median(plot_data[i]) / gtex_median if gtex_median and gtex_median > 0 else float('inf')
            fold_changes.append(f"FC={fc:.1f}")
        
        if fold_changes:
            ax1.text(0.98, 0.95, f"Normal Brain median TPM: {gtex_median:.2f}\n" + "\n".join(fold_changes), 
                    transform=ax1.transAxes, ha='right', va='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # B: 密度图
    ax2 = fig.add_subplot(gs[1, :2])
    for i, (data, label, c) in enumerate(zip(plot_data, labels, colors)):
        if len(data) > 1:
            ax2.hist(np.log2(data + 0.01), bins=30, alpha=0.5, color=c, label=label, density=True)
    ax2.set_xlabel("log2(TPM + 0.01)", fontsize=12)
    ax2.set_ylabel("Density", fontsize=12)
    ax2.set_title("B. ZP3 Expression Distribution (log2 scale)", fontsize=13, fontweight='bold')
    ax2.legend()
    
    # C: 统计摘要表
    ax3 = fig.add_subplot(gs[1, 2])
    ax3.axis('off')
    summary_lines = ["Group      n    Median  Mean"]
    summary_lines.append("-" * 35)
    for data, label in zip(plot_data, labels):
        lab_short = label.split("\n")[0]
        summary_lines.append(f"{lab_short:12s} {len(data):4d}  {np.median(data):6.2f}  {np.mean(data):6.2f}")
    ax3.text(0.05, 0.95, "\n".join(summary_lines), transform=ax3.transAxes, fontsize=10,
            va='top', fontfamily='monospace')
    ax3.set_title("C. Expression Statistics", fontsize=12, fontweight='bold')
    
    # D: 统计检验
    ax4 = fig.add_subplot(gs[2, :])
    tests = []
    ax4.axis('off')
    
    test_text = "Statistical Tests (Mann-Whitney U)\n" + "=" * 50 + "\n"
    
    # GBM vs GTEx Brain
    if len(plot_data) >= 2 and len(plot_data[0]) > 0 and len(plot_data[1]) > 0:
        # GTEx 是 plot_data[0], GBM 是 plot_data[1]
        stat, p = stats.mannwhitneyu(plot_data[1], plot_data[0], alternative='two-sided')
        fc = np.median(plot_data[1]) / (np.median(plot_data[0]) + 1e-10)
        test_text += f"GBM vs GTEx Brain:  FC={fc:.1f}, p={p:.2e}  {'***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'}\n"
    
    if len(plot_data) >= 2 and len(plot_data[0]) > 0 and len(plot_data[-1]) > 0:
        stat, p = stats.mannwhitneyu(plot_data[-1], plot_data[0], alternative='two-sided')
        fc = np.median(plot_data[-1]) / (np.median(plot_data[0]) + 1e-10)
        test_text += f"LGG vs GTEx Brain:  FC={fc:.1f}, p={p:.2e}  {'***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'}\n"
    
    test_text += f"\nGTEx Brain median: {np.median(plot_data[0]):.4f}"
    if len(plot_data) >= 2:
        test_text += f"\nGBM median: {np.median(plot_data[1]):.4f} (FC={np.median(plot_data[1])/(np.median(plot_data[0])+1e-10):.1f}x)"
    if len(plot_data) >= 3:
        test_text += f"\nLGG median: {np.median(plot_data[2]):.4f} (FC={np.median(plot_data[2])/(np.median(plot_data[0])+1e-10):.1f}x)"
    
    ax4.text(0.05, 0.95, test_text, transform=ax4.transAxes, fontsize=11,
            va='top', fontfamily='monospace')
    ax4.set_title("D. Statistical Comparison", fontsize=12, fontweight='bold')
    
    # 保存
    fig_path = os.path.join(OUTPUT_DIR, "fig_zp3_normal_vs_tumor.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"  图表: {fig_path}")
    plt.close()
    
    return fig_path

def save_results(gtex_brain, df_tcga, fold_changes):
    """保存结果 CSV"""
    # GTEx 脑数据
    if gtex_brain is not None:
        gtex_brain.to_csv(os.path.join(OUTPUT_DIR, "gtex_brain_zp3_expression.csv"))
    
    # TCGA
    if df_tcga is not None:
        df_tcga.to_csv(os.path.join(OUTPUT_DIR, "tcga_zp3_expression.csv"))
    
    # 摘要
    summary_path = os.path.join(OUTPUT_DIR, "normal_vs_tumor_summary.csv")
    # 收集统计
    rows = []
    for grp, data in [("GTEx_Brain", gtex_brain), ("TCGA_GBM", df_tcga[df_tcga["Project"]=="GBM"] if df_tcga is not None else None), ("TCGA_LGG", df_tcga[df_tcga["Project"]=="LGG"] if df_tcga is not None else None)]:
        if data is not None and len(data) > 0:
            col = "ZP3_TPM" if "ZP3_TPM" in data.columns else data.columns[0]
            vals = data[col].dropna()
            rows.append({"Group": grp, "n": len(vals), "Mean": vals.mean(), "Median": vals.median(), "Std": vals.std()})
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    print(f"  摘要: {summary_path}")

def main():
    print("=" * 60)
    print("GTEx 正常脑对照: ZP3 正常 vs 肿瘤表达差异")
    print("=" * 60)
    
    # 1. GTEx phenotype
    phenotype_path = download_gtex_phenotype()
    brain_samples, brain_col = load_phenotype(phenotype_path)
    
    # 2. 提取 ZP3 表达
    import zlib  # for gzip decompression in download
    df_brain = extract_gtex_zp3_expression("zp3_toil_isoform_tpm.tsv", brain_samples)
    
    # 3. TCGA
    df_tcga = load_tcga_zp3_expression()
    
    # 4. 对比
    fc = None
    if df_brain is not None and df_tcga is not None:
        gtex_med = df_brain["ZP3_TPM"].median() if "ZP3_TPM" in df_brain.columns else 0
        gbm_med = df_tcga[df_tcga["Project"]=="GBM"]["ZP3_TPM"].median() if len(df_tcga[df_tcga["Project"]=="GBM"])>0 else 0
        fc = {"GBM_vs_GTEx": gbm_med/gtex_med if gtex_med>0 else float('inf')}
    
    fig_path = create_comparison(df_brain, df_tcga, None)
    save_results(df_brain, df_tcga, fc)
    
    print(f"\n✓ GTEx 正常脑对照分析完成!")
    print(f"  产物目录: {OUTPUT_DIR}")

if __name__ == "__main__":
    import zlib  # early import for download function
    main()
