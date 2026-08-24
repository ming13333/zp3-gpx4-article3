#!/usr/bin/env python3
"""
GTEx normal brain control analysis: ZP3 normal vs tumor expression difference

Data sources:
1. GTEX_phenotype (Xena Toil hub) — GTEx sample tissue mapping
2. Toil isoform TPM (downloading) — ZP3 gene total TPM = sum(isoform TPM)
3. cBioPortal GBM/LGG gene expression (already available)
4. HPA protein expression (already available)

Output: ZP3 normal brain vs GBM vs LGG expression comparison + fold change
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
# 1. Download GTEX phenotype
# ======================
def download_gtex_phenotype():
    """Download GTEx phenotype data (sample -> tissue)"""
    print("[1/5] Downloading GTEx phenotype data...")
    url = "https://toil.xenahubs.net/download/GTEX_phenotype.gz"
    local_path = os.path.join(DIR, "GTEX_phenotype.tsv")
    
    if os.path.exists(local_path) and os.path.getsize(local_path) > 1000:
        print(f"  Already exists: {local_path}")
        return local_path
    
    print("  Streaming download + decompression...")
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
    
    # Write to file
    with open(local_path, "w", encoding="utf-8") as f:
        for line in text_lines:
            f.write(line + "\n")
    
    print(f"  Saved: {local_path} ({len(text_lines)} lines)")
    return local_path

def load_phenotype(path):
    """Load GTEx phenotype, filter brain tissue samples"""
    print("\n[2/5] Filter normal brain tissue samples...")
    df = pd.read_csv(path, sep="\t", index_col=0)
    print(f"  Total GTEx samples: {len(df)}")
    print(f"  Columns: {list(df.columns)}")
    
    # Find tissue-related columns
    tissue_cols = [c for c in df.columns if "tissue" in c.lower() or "site" in c.lower() or "body" in c.lower()]
    if tissue_cols:
        print(f"  Tissue-related columns: {tissue_cols[:5]}")
    
    # Find brain tissue
    brain_col = None
    for col in df.columns:
        # Check if this column contains brain
        if df[col].dtype == object:
            brain_hits = df[col].astype(str).str.lower().str.contains("brain|brain_|cerebr|hippocamp|cortex|frontal|temporal|cerebell", na=False)
            if brain_hits.sum() > 10:
                brain_col = col
                break
    
    if brain_col:
        brain_terms = df[brain_col].astype(str).str.lower()
        is_brain = brain_terms.str.contains("brain|brain_|cerebr|hippocamp|cortex|frontal|temporal|cerebell", na=False)
        brain_samples = df[is_brain]
        print(f"\n  Brain tissue column: {brain_col}")
        print(f"  Brain tissue samples: {len(brain_samples)}")
        
        # Subdivide brain regions
        if df[brain_col].dtype == object:
            brain_regions = brain_samples[brain_col].value_counts().head(15)
            print(f"  Brain region distribution:\n{brain_regions.to_string()}")
        
        return brain_samples, brain_col
    else:
        # Try all columns
        print("  No dedicated tissue column found, printing first few sample rows...")
        print(f"  First 5 rows: {df.head(3)}")
        return df, None

# ======================
# 2. Extract GTEx sample expression from downloaded ZP3 isoform data
# ======================
def extract_gtex_zp3_expression(isoform_tsv, brain_samples):
    """Extract total ZP3 expression of GTEx brain samples from Toil isoform TPM"""
    print("\n[3/5] Extracting ZP3 expression from GTEx brain samples...")
    
    isoform_path = os.path.join(DIR, isoform_tsv)
    if not os.path.exists(isoform_path):
        print(f"  ⚠ isoform data not downloaded yet: {isoform_path}")
        return None
    
    df_iso = pd.read_csv(isoform_path, sep="\t", index_col=0)
    df_iso = df_iso.T  # samples × transcripts
    
    # Filter GTEx samples
    gtex_idx = [s for s in df_iso.index if str(s).startswith("GTEX")]
    df_gtex = df_iso.loc[gtex_idx]
    
    # Calculate total TPM (log2 scale → TPM)
    for col in df_gtex.columns:
        df_gtex[f"{col}_tpm"] = 2**df_gtex[col] - 0.001
    
    tpm_cols = [c for c in df_gtex.columns if c.endswith("_tpm")]
    df_gtex["ZP3_TPM"] = df_gtex[tpm_cols].sum(axis=1)
    df_gtex["ZP3_TPM_log2"] = np.log2(df_gtex["ZP3_TPM"] + 0.001)
    
    # Match brain tissue samples
    gtex_ids = set(brain_samples.index) if brain_samples is not None else set()
    brain_in_gtex = df_gtex.index.intersection(gtex_ids)
    print(f"   Brain tissue samples in isoform data: {len(brain_in_gtex)}/{len(gtex_ids)}")
    
    df_brain = df_gtex.loc[brain_in_gtex] if len(brain_in_gtex) > 0 else df_gtex
    
    return df_brain, df_gtex

# ======================
# 3. Load cBioPortal GBM/LGG ZP3 expression
# ======================
def load_tcga_zp3_expression():
    """Load TCGA GBM/LGG ZP3 from existing cBioPortal data"""
    print("\n[4/5] Loading TCGA GBM/LGG ZP3 expression (cBioPortal)...")
    
    # Read from existing file or create new
    # Use cBioPortal API to obtain
    results = {}
    
    # ZP3 Entrez = 7784 (audited 2026-08-24; was erroneously 8277 = SP5)
    for study, entrez_id in [("gbm_tcga", 7784), ("lgg_tcga", 7784)]:
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
                print(f"  {study}: {len(values)} samples")
        except Exception as e:
            print(f"  {study}: failed {e}")
    
    # Merge
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
    
    print(f"  TCGA total: {len(df_tcga)} samples")
    return df_tcga

# ======================
# 4. Integration + visualization
# ======================
def create_comparison(df_brain, df_tcga, gtex_all):
    """Create normal vs tumor comparison visualization + statistics"""
    print("\n[5/5] Normal vs tumor comparison analysis...")
    
    fig = plt.figure(figsize=(14, 12))
    gs = plt.matplotlib.gridspec.GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # A: Boxplot — GTEx Brain vs GBM vs LGG
    ax1 = fig.add_subplot(gs[0, :])
    
    plot_data = []
    labels = []
    colors = []
    
    # GTEx brain
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
    
    # Boxplot
    bp = ax1.boxplot(plot_data, patch_artist=True, widths=0.5, showfliers=True, flierprops={'alpha': 0.3, 'markersize': 3})
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    
    ax1.set_xticklabels(labels, fontsize=12)
    ax1.set_ylabel("ZP3 Expression (TPM/FPM)", fontsize=13)
    ax1.set_title("A. ZP3 Expression: Normal Brain vs GBM vs LGG", fontsize=14, fontweight='bold')
    ax1.set_yscale("log")
    
    # Calculate fold change
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
    
    # B: Density plot
    ax2 = fig.add_subplot(gs[1, :2])
    for i, (data, label, c) in enumerate(zip(plot_data, labels, colors)):
        if len(data) > 1:
            ax2.hist(np.log2(data + 0.01), bins=30, alpha=0.5, color=c, label=label, density=True)
    ax2.set_xlabel("log2(TPM + 0.01)", fontsize=12)
    ax2.set_ylabel("Density", fontsize=12)
    ax2.set_title("B. ZP3 Expression Distribution (log2 scale)", fontsize=13, fontweight='bold')
    ax2.legend()
    
    # C: Statistical summary table
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
    
    # D: Statistical test
    ax4 = fig.add_subplot(gs[2, :])
    tests = []
    ax4.axis('off')
    
    test_text = "Statistical Tests (Mann-Whitney U)\n" + "=" * 50 + "\n"
    
    # GBM vs GTEx Brain
    if len(plot_data) >= 2 and len(plot_data[0]) > 0 and len(plot_data[1]) > 0:
        # GTEx is plot_data[0], GBM is plot_data[1]
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
    
    # Save
    fig_path = os.path.join(OUTPUT_DIR, "fig_zp3_normal_vs_tumor.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"  Figure: {fig_path}")
    plt.close()
    
    return fig_path

def save_results(gtex_brain, df_tcga, fold_changes):
    """Save results CSV"""
    # GTEx brain data
    if gtex_brain is not None:
        gtex_brain.to_csv(os.path.join(OUTPUT_DIR, "gtex_brain_zp3_expression.csv"))
    
    # TCGA
    if df_tcga is not None:
        df_tcga.to_csv(os.path.join(OUTPUT_DIR, "tcga_zp3_expression.csv"))
    
    # Summary
    summary_path = os.path.join(OUTPUT_DIR, "normal_vs_tumor_summary.csv")
    # Collect statistics
    rows = []
    for grp, data in [("GTEx_Brain", gtex_brain), ("TCGA_GBM", df_tcga[df_tcga["Project"]=="GBM"] if df_tcga is not None else None), ("TCGA_LGG", df_tcga[df_tcga["Project"]=="LGG"] if df_tcga is not None else None)]:
        if data is not None and len(data) > 0:
            col = "ZP3_TPM" if "ZP3_TPM" in data.columns else data.columns[0]
            vals = data[col].dropna()
            rows.append({"Group": grp, "n": len(vals), "Mean": vals.mean(), "Median": vals.median(), "Std": vals.std()})
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    print(f"  Summary: {summary_path}")

def main():
    print("=" * 60)
    print("GTEx normal brain control: ZP3 normal vs tumor expression differences")
    print("=" * 60)
    
    # 1. GTEx phenotype
    phenotype_path = download_gtex_phenotype()
    brain_samples, brain_col = load_phenotype(phenotype_path)
    
    # 2. Extract ZP3 expression
    import zlib  # for gzip decompression in download
    df_brain = extract_gtex_zp3_expression("zp3_toil_isoform_tpm.tsv", brain_samples)
    
    # 3. TCGA
    df_tcga = load_tcga_zp3_expression()
    
    # 4. Comparison
    fc = None
    if df_brain is not None and df_tcga is not None:
        gtex_med = df_brain["ZP3_TPM"].median() if "ZP3_TPM" in df_brain.columns else 0
        gbm_med = df_tcga[df_tcga["Project"]=="GBM"]["ZP3_TPM"].median() if len(df_tcga[df_tcga["Project"]=="GBM"])>0 else 0
        fc = {"GBM_vs_GTEx": gbm_med/gtex_med if gtex_med>0 else float('inf')}
    
    fig_path = create_comparison(df_brain, df_tcga, None)
    save_results(df_brain, df_tcga, fc)
    
    print(f"\n✓ GTEx normal brain control analysis completed!")
    print(f"  Output directory: {OUTPUT_DIR}")

if __name__ == "__main__":
    import zlib  # early import for download function
    main()
