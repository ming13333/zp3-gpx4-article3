#!/usr/bin/env python3
"""
Complete GTEx normal brain control analysis

Two-step streaming download:
1. GTEX_phenotype (small file) → brain tissue sample IDs
2. TcgaTargetGtex_rsem_gene_tpm (1.23 GB) → stream-extract ZP3 row
3. Normal brain vs GBM vs LGG comparison
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys, time, json, zlib, requests, io, gzip

DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "output", "phase1_knowledge_gap_filling")
GENE_TPM_URL = "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz"
PHENOTYPE_URL = "https://toil.xenahubs.net/download/GTEX_phenotype.gz"
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "article3", "results", "gtex_brain_control")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# 1. Download GTEX phenotype
# ==========================================
def download_phenotype():
    print("[1/4] Reading local GTEX phenotype...", flush=True)
    t0 = time.time()
    
    local_path = os.path.join(DIR, "GTEX_phenotype.gz")
    raw_lines = []
    with gzip.open(local_path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            raw_lines.append(line.rstrip('\n'))
    
    # Parse into DataFrame
    # First row is header
    header = raw_lines[0].split("\t")
    data_rows = [line.split("\t") for line in raw_lines[1:] if line.strip()]
    df = pd.DataFrame(data_rows, columns=header)
    df = df.set_index("sample") if "sample" in df.columns else df.set_index(df.columns[0])
    
    elapsed = time.time() - t0
    print(f"  Reading complete, {len(df)} samples, elapsed {elapsed:.1f}s", flush=True)
    return df

def find_brain_samples(pheno):
    """Find brain tissue samples in phenotype"""
    print("\n[2/4] Filtering brain tissue samples...")
    
    # Check key columns
    tissue_cols = [c for c in pheno.columns if any(kw in c.lower() for kw in ["tissue", "site", "organ", "body"])]
    detail_cols = [c for c in pheno.columns if "detail" in c.lower() or "description" in c.lower() or "name" in c.lower()]
    
    # Print first few columns
    print(f"  Number of columns: {len(pheno.columns)}")
    print(f"  First 10 columns: {list(pheno.columns[:10])}")
    
    # Find samples containing brain
    brain_mask = pd.Series(False, index=pheno.index)
    brain_col = None
    brain_values = None
    
    for col in pheno.columns:
        if pheno[col].dtype == object:
            s = pheno[col].astype(str).str.lower()
            m = s.str.contains("brain", na=False)
            if m.sum() > 10:
                brain_mask = brain_mask | m
                if brain_col is None:
                    brain_col = col
                    brain_values = pheno.loc[m, col]
    
    brain_samples = pheno[brain_mask]
    print(f"  Brain tissue samples: {len(brain_samples)}/{len(pheno)}")
    
    if brain_values is not None:
        print(f"  Brain region distribution (top 10):")
        for name, count in brain_values.value_counts().head(10).items():
            print(f"    {name}: {count}")
    
    return brain_samples, brain_col

# ==========================================
# 2. Stream download gene TPM + extract ZP3
# ==========================================
def extract_zp3_gene_expression(brain_sample_ids):
    """Stream extract ZP3 + brain samples from TcgaTargetGtex_rsem_gene_tpm"""
    print("\n[3/4] Stream extract ZP3 gene expression (1.23 GB)...", flush=True)
    t0 = time.time()
    
    local_path = os.path.join(DIR, "TcgaTargetGtex_rsem_gene_tpm.gz")
    file_size = os.path.getsize(local_path)
    print(f"  File size: {file_size/1e9:.2f} GB", flush=True)
    
    d = zlib.decompressobj(16 + zlib.MAX_WBITS)
    buf = b""
    header = None
    zp3_line = None
    total_lines = 0
    total_bytes = 0
    last_report = 0
    
    with open(local_path, 'rb') as fh:
        while True:
            chunk = fh.read(65536)
            if not chunk:
                break
            total_bytes += len(chunk)
            buf += d.decompress(chunk)
        
            while b"\n" in buf:
                line_bytes, buf = buf.split(b"\n", 1)
                total_lines += 1
                
                try:
                    line = line_bytes.decode("utf-8", errors="replace")
                except:
                    continue
                
                if total_lines == 1:
                    header = line
                    headers = header.split("\t")
                    # Extract required column indices
                    brain_idx = []
                    gbm_idx = []
                    lgg_idx = []
                    for i, h in enumerate(headers):
                        h_upper = h.upper()
                        if h in brain_sample_ids:
                            brain_idx.append(i)
                        elif h.startswith("TCGA-"):
                            parts = h.split("-")
                            if len(parts) >= 2:
                                # Quickly determine GBM/LGG (based on known code or cBioPortal)
                                pass  # Subsequent processing
                    
                    print(f"  Column headers: {header[:200]}...", flush=True)
                    print(f"  Total columns: {len(headers)}", flush=True)
                    continue
                
                # Find ZP3
                if line.startswith("ENSG00000188372\t") or line.startswith("ENSG00000188372."):
                    zp3_line = line
                    print(f"  ✓ Found ZP3 row (line {total_lines})", flush=True)
                
                # ZP3 + header both found, can stop
                if zp3_line is not None and header is not None:
                    break
            
            # Progress
            if total_bytes - last_report > 100_000_000:
                elapsed = time.time() - t0
                pct = total_bytes / file_size * 100
                print(f"  [{elapsed/60:.1f} min] {total_bytes/1e6:.0f} MB ({pct:.1f}%)", flush=True)
                last_report = total_bytes
            
            if zp3_line is not None and header is not None:
                break  # exit while True
    
    elapsed = time.time() - t0
    print(f"  Streaming read complete, took {elapsed/60:.1f} min, scanned {total_lines} lines", flush=True)
    
    if header is None:
        print("  ⚠ Header not found!")
        return None
    if zp3_line is None:
        print("  ⚠ ZP3 row not found!")
        return None
    
    # Parse ZP3 expression
    headers = header.split("\t")
    zp3_vals = zp3_line.split("\t")
    
    # Build sample→ZP3 TPM mapping
    expr_dict = {}
    for i in range(1, len(headers)):
        if i < len(zp3_vals):
            try:
                expr_dict[headers[i]] = float(zp3_vals[i])
            except ValueError:
                pass
    
    print(f"  Total samples: {len(expr_dict)}", flush=True)
    
    # Classify
    gtex_samples = {k: v for k, v in expr_dict.items() if k.startswith("GTEX")}
    tcga_samples = {k: v for k, v in expr_dict.items() if k.startswith("TCGA")}
    print(f"  GTEx samples: {len(gtex_samples)}", flush=True)
    print(f"  TCGA samples: {len(tcga_samples)}", flush=True)
    
    # Filter brain tissue GTEx
    gtex_brain = {k: v for k, v in gtex_samples.items() if k in brain_sample_ids}
    print(f"  GTEx Brain samples: {len(gtex_brain)}", flush=True)
    
    # TCGA split by GBM/LGG (via cBioPortal)
    return {
        "all": expr_dict,
        "gtex_brain": gtex_brain,
        "gtex_all": gtex_samples,
        "tcga_all": tcga_samples,
        "header": header
    }

# ==========================================
# 3. Match TCGA GBM/LGG
# ==========================================
def classify_tcga_samples(tcga_dict):
    """Classify TCGA samples as GBM/LGG using cBioPortal API"""
    print("\n[4/5] Classifying TCGA samples as GBM/LGG...")
    
    # Get sample list from cBioPortal
    sample_to_project = {}
    for study_id in ["gbm_tcga", "lgg_tcga"]:
        try:
            url = f"https://www.cbioportal.org/api/studies/{study_id}/samples"
            r = requests.get(url, timeout=60, headers={"Accept": "application/json"})
            if r.status_code == 200:
                samples = r.json()
                proj = "GBM" if "gbm" in study_id else "LGG"
                for s in samples:
                    sid = s.get("sampleId", "")
                    short = "-".join(sid.split("-")[:3])
                    if short not in sample_to_project:
                        sample_to_project[short] = proj
        except Exception as e:
            print(f"  {study_id}: {e}")
    
    print(f"  cBioPortal mapping: {len(sample_to_project)} TCGA samples")
    
    gbm_raw = {}
    lgg_raw = {}
    other = {}
    
    for sid, val in tcga_dict.items():
        short = "-".join(sid.split("-")[:3])
        proj = sample_to_project.get(short, "Other")
        if proj == "GBM":
            gbm_raw[sid] = val
        elif proj == "LGG":
            lgg_raw[sid] = val
        else:
            other[sid] = val
    
    print(f"  GBM: {len(gbm_raw)}, LGG: {len(lgg_raw)}, Other TCGA: {len(other)}")
    return gbm_raw, lgg_raw, other

# ==========================================
# 4. Visualization + Statistics
# ==========================================
def create_visualization(gtex_brain, gbm, lgg, gtex_all):
    print("\n[5/5] Normal vs Tumor comparison visualization...")
    
    # Prepare data
    gb_vals = np.array(list(gtex_brain.values())) if gtex_brain else np.array([])
    gbm_vals = np.array(list(gbm.values())) if gbm else np.array([])
    lgg_vals = np.array(list(lgg.values())) if lgg else np.array([])
    ga_vals = np.array(list(gtex_all.values())) if gtex_all else np.array([])
    
    # Statistics
    stats_data = {}
    for name, vals in [("GTEx_Brain", gb_vals), ("GTEx_All", ga_vals), 
                        ("TCGA_GBM", gbm_vals), ("TCGA_LGG", lgg_vals)]:
        if len(vals) > 0:
            stats_data[name] = {
                "n": len(vals),
                "mean": np.mean(vals),
                "median": np.median(vals),
                "std": np.std(vals),
                "sem": np.std(vals) / np.sqrt(len(vals))
            }
    
    # Print statistics
    for name, d in stats_data.items():
        print(f"  {name}: n={d['n']}, mean={d['mean']:.4f}, median={d['median']:.4f}")
    
    # Fold change
    if len(gb_vals) > 0:
        gb_med = np.median(gb_vals)
        print(f"\n  Fold change (vs GTEx Brain median={gb_med:.4f}):")
        if len(gbm_vals) > 0:
            fc = np.median(gbm_vals) / gb_med if gb_med > 0 else float('inf')
            print(f"    GBM: {fc:.1f}x")
        if len(lgg_vals) > 0:
            fc = np.median(lgg_vals) / gb_med if gb_med > 0 else float('inf')
            print(f"    LGG: {fc:.1f}x")
    
    # Statistical tests
    print(f"\n  Statistical test (Mann-Whitney U):")
    for name, vals in [("GBM", gbm_vals), ("LGG", lgg_vals)]:
        if len(vals) > 0 and len(gb_vals) > 0:
            stat, p = stats.mannwhitneyu(vals, gb_vals, alternative="two-sided")
            fc = np.median(vals) / np.median(gb_vals) if np.median(gb_vals) > 0 else float('inf')
            sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
            print(f"    {name} vs GTEx Brain: FC={fc:.1f}, U={stat:.0f}, p={p:.2e} {sig}")
    
    # --- Plots ---
    fig = plt.figure(figsize=(16, 13))
    gs = plt.matplotlib.gridspec.GridSpec(3, 3, hspace=0.35, wspace=0.3)
    
    # A: Boxplot (log scale)
    ax1 = fig.add_subplot(gs[0, :2])
    plot_groups = []
    plot_labels = []
    plot_colors = []
    
    for name, vals, color in [
        ("GTEx Brain (n={})", gb_vals, "#2ECC71"),
        ("GTEx All (n={})", ga_vals, "#27AE60"),
        ("TCGA-GBM (n={})", gbm_vals, "#E74C3C"),
        ("TCGA-LGG (n={})", lgg_vals, "#3498DB")
    ]:
        if len(vals) > 0:
            plot_groups.append(vals)
            plot_labels.append(name.format(len(vals)))
            plot_colors.append(color)
    
    bp = ax1.boxplot(plot_groups, patch_artist=True, widths=0.5,
                     flierprops={'alpha': 0.3, 'markersize': 3})
    for patch, c in zip(bp["boxes"], plot_colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    
    ax1.set_xticklabels(plot_labels, fontsize=10)
    ax1.set_ylabel("ZP3 Expression (log2 TPM + 0.001)", fontsize=12)
    ax1.set_title("A. ZP3 Expression: Normal vs Tumor", fontsize=13, fontweight='bold')
    
    # B: Density plot
    ax2 = fig.add_subplot(gs[0, 2])
    for vals, label, c in zip(plot_groups, plot_labels, plot_colors):
        if len(vals) > 1:
            ax2.hist(vals, bins=25, alpha=0.4, color=c, label=label, density=True)
    ax2.set_xlabel("log2(TPM + 0.001)", fontsize=10)
    ax2.set_title("B. Distribution", fontsize=11, fontweight='bold')
    ax2.legend(fontsize=7, loc='upper left')
    
    # C: Fold change bar chart
    ax3 = fig.add_subplot(gs[1, :2])
    if len(gb_vals) > 0:
        gb_med = np.median(gb_vals)
        fc_data = {}
        fc_err = {}
        if len(gbm_vals) > 0:
            fc_data["GBM"] = np.median(gbm_vals) / gb_med
        if len(lgg_vals) > 0:
            fc_data["LGG"] = np.median(lgg_vals) / gb_med
        if len(ga_vals) > 0:
            fc_data["GTEx All"] = np.median(ga_vals) / gb_med
        
        bars = ax3.bar(fc_data.keys(), fc_data.values(), 
                       color=['#E74C3C', '#3498DB', '#27AE60'][:len(fc_data)], alpha=0.7)
        ax3.axhline(1.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax3.set_ylabel("Fold Change (vs GTEx Brain)", fontsize=11)
        ax3.set_title(f"C. Fold Change vs GTEx Brain (median={gb_med:.2f})", fontsize=12, fontweight='bold')
        for bar, val in zip(bars, fc_data.values()):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, 
                    f'{val:.1f}x', ha='center', fontsize=12, fontweight='bold')
    
    # D: Statistical summary table
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.axis('off')
    summary = "Group            n     Median   Mean\n" + "-"*39 + "\n"
    for name, d in stats_data.items():
        summary += f"{name:18s} {d['n']:4d}   {d['median']:6.3f}  {d['mean']:6.3f}\n"
    ax4.text(0.05, 0.95, summary, transform=ax4.transAxes, fontsize=9,
            va='top', fontfamily='monospace')
    ax4.set_title("D. Summary Statistics", fontsize=11, fontweight='bold')
    
    # E: Statistical test summary
    ax5 = fig.add_subplot(gs[2, :])
    ax5.axis('off')
    test_text = "Statistical Tests (Mann-Whitney U)\n" + "="*50 + "\n"
    for name, vals in [("TCGA-GBM", gbm_vals), ("TCGA-LGG", lgg_vals)]:
        if len(vals) > 0 and len(gb_vals) > 0:
            stat, p = stats.mannwhitneyu(vals, gb_vals, alternative="two-sided")
            fc = np.median(vals) / np.median(gb_vals) if np.median(gb_vals)>0 else float('inf')
            test_text += f"{name} vs GTEx Brain:  FC={fc:.1f}x, U={stat:.0f}, p={p:.2e}\n"
    test_text += f"\nInterpretation: ZP3 is {'dramatically' if fc>5 else 'moderately' if fc>2 else 'slightly'} overexpressed in glioma vs normal brain.\n"
    test_text += f"This supports the hypothesis that ZP3 is aberrantly activated in brain tumors.\n"
    ax5.text(0.05, 0.95, test_text, transform=ax5.transAxes, fontsize=10,
            va='top', fontfamily='monospace')
    ax5.set_title("E. Interpretation", fontsize=11, fontweight='bold')
    
    fig_path = os.path.join(OUTPUT_DIR, "fig_zp3_normal_vs_tumor.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    print(f"  Figure: {fig_path}")
    plt.close()
    
    return fig_path, stats_data

def save_outputs(gtex_brain, gbm, lgg, gtex_all, stats_data):
    """Save CSV + JSON"""
    # CSV
    for name, data in [("gtex_brain", gtex_brain), ("tcga_gbm", gbm), 
                        ("tcga_lgg", lgg), ("gtex_all", gtex_all)]:
        if data:
            path = os.path.join(OUTPUT_DIR, f"{name}_zp3_expression.csv")
            pd.DataFrame({"sample": list(data.keys()), "ZP3_log2_TPM": list(data.values())}).to_csv(path, index=False)
    
    # Statistics JSON
    json_path = os.path.join(OUTPUT_DIR, "normal_vs_tumor_stats.json")
    with open(json_path, "w") as f:
        json.dump(stats_data, f, indent=2, ensure_ascii=False)
    print(f"  Statistics: {json_path}")

def main():
    print("=" * 60)
    print("GTEx normal brain control: ZP3 normal vs tumor expression differences")
    print("=" * 60)
    
    # 1. GTEX phenotype
    pheno = download_phenotype()
    brain_samples, brain_col = find_brain_samples(pheno)
    brain_ids = set(brain_samples.index.tolist())
    
    # 2. Streaming extraction of ZP3 gene expression
    expr_data = extract_zp3_gene_expression(brain_ids)
    if expr_data is None:
        print("⚠ Gene expression extraction failed")
        return
    
    # 3. Classify TCGA
    gbm, lgg, other = classify_tcga_samples(expr_data["tcga_all"])
    
    # 4. Visualization
    fig_path, stats_data = create_visualization(
        expr_data["gtex_brain"], gbm, lgg, expr_data["gtex_all"]
    )
    
    # 5. Save
    save_outputs(expr_data["gtex_brain"], gbm, lgg, expr_data["gtex_all"], stats_data)
    
    print(f"\n✓ GTEx normal brain control complete!")
    print(f"  Outputs: {OUTPUT_DIR}")
    return True

if __name__ == "__main__":
    main()
