# -*- coding: utf-8 -*-
"""
Article 3 Supplementary ② — GSEA functional annotation of FL-high vs RI-high
=====================================================
Purpose: Examine functional implications of isoform switching — samples with high FL full-length proportion vs RI (retained-intron)
samples with high proportion, which pathways are enriched in differentially expressed genes (focus: immunosuppressive pathways).

Design:
  1. Samples: GBM+LGG combined glioma (TCGA), grouped by PSI
     - FL_high: FL-PSI ≥ upper tertile and RI-PSI ≤ lower tertile (pure FL group)
     - RI_high: RI-PSI ≥ upper tertile and FL-PSI ≤ lower tertile (pure RI group)
  2. Differential expression: per-gene Wilcoxon rank-sum test between groups → z-score as ranking metric
     (z = normal approximation statistic, positive = upregulated in FL_high)
  3. GSEA: gseapy.prerank, gene sets MSigDB Hallmark_2020 + C7 immune
  4. Output: enrichment table + enrichment plots for key immunosuppressive pathways

Outputs: zp3_gsea_results/
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

FL = "ENST00000336517.8"   # canonical full-length
RI = "ENST00000466960.5"   # retained-intron
HALLMARK = "MSigDB_Hallmark_2020"
C7_IMMUNE = "MSigDB_Immunologic_Signatures"   # skip if unavailable
ENSG_SYMBOL_CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "phase1_knowledge_gap_filling", "ensg_symbol_map.json")


def read_psi():
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    return psi


def load_expression(target_ensg_set, sample_subset):
    """Stream-read target genes × target samples from TPM (returns genes × samples)."""
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
                print(f"    scanned {n} lines...")
    return pd.DataFrame(rows, index=sample_list).T  # genes × samples


def main():
    print("=== Article 3 Supplement ②: FL-high vs RI-high GSEA ===\n")

    # 1. PSI + cancer type
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
    print(f"Glioma samples (GBM+LGG): {len(gl)}")

    sub = psi.loc[gl, [FL, RI]].dropna()
    print(f"Samples with PSI values: {len(sub)}")

    # 2. Tertile grouping (pure groups)
    fl_q = sub[FL].quantile(0.67)
    ri_q = sub[RI].quantile(0.67)
    fl_low_q = sub[FL].quantile(0.33)
    ri_low_q = sub[RI].quantile(0.33)
    print(f"FL tertiles: {sub[FL].quantile(0.33):.3f} / {fl_q:.3f}")
    print(f"RI tertiles: {sub[RI].quantile(0.33):.3f} / {ri_q:.3f}")

    fl_high = sub[(sub[FL] >= fl_q) & (sub[RI] <= ri_low_q)].index
    ri_high = sub[(sub[RI] >= ri_q) & (sub[FL] <= fl_low_q)].index
    print(f"FL_high (pure FL): {len(fl_high)} | RI_high (pure RI): {len(ri_high)}")
    if len(fl_high) < 15 or len(ri_high) < 15:
        print("!! Insufficient group samples, relaxing criteria")
        fl_high = sub[sub[FL] >= fl_q].index
        ri_high = sub[sub[RI] >= ri_q].index
        print(f"  After relaxing: FL_high={len(fl_high)} RI_high={len(ri_high)}")

    # 3. Read full expression (both groups)
    print("Reading TPM (both groups)...")
    all_samp = sorted(set(fl_high) | set(ri_high))
    # Target gene set = immune-related (speed up reading; expand when GSEA needs full gene space)
    # Read all genes first (~60000 rows, streaming feasible)
    ensg_target = None  # None = read all genes
    # But reading all genes uses large memory; compromise: skipping low-expression gene rows lowers quality, so read all
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
                print(f"    Scanned {n} lines...")
    expr = pd.DataFrame(rows, index=all_samp).T  # genes × samples
    print(f"Expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

    # 4. Differential expression (Wilcoxon, two groups)
    print("Differential expression test (Wilcoxon)...")
    g_fl = expr[fl_high]
    g_ri = expr[ri_high]
    genes = list(expr.index)
    recs = []
    for g in genes:
        a = g_fl.loc[g].values.astype(float)
        b = g_ri.loc[g].values.astype(float)
        # TPM file is log2(TPM+ε); filter extremely low expression (log2 mean < -3 ≈ TPM<0.12)
        if np.nanmean(a) < -3.0 and np.nanmean(b) < -3.0:
            continue  # filter low expression
        try:
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            # z approximation
            na, nb = len(a), len(b)
            mu = na * nb / 2
            sd = np.sqrt(na * nb * (na + nb + 1) / 12)
            z = (u - mu) / sd if sd > 0 else 0.0
            recs.append({"gene": g, "z": z, "p": p,
                         "mean_FL": float(np.nanmean(a)), "mean_RI": float(np.nanmean(b))})
        except Exception:
            continue
    de = pd.DataFrame(recs)
    print(f"DE gene count: {len(de)}")

    # Gene name: ensg -> symbol. Prefer local whole-genome mapping cache; otherwise fall back to ensg_map.
    ensg2sym = {}
    if os.path.exists(ENSG_SYMBOL_CACHE):
        with open(ENSG_SYMBOL_CACHE) as f:
            ensg2sym = json.load(f)
        print(f"Loaded whole-genome symbol mapping: {len(ensg2sym)} entries")
    else:
        print("!! No whole-genome mapping cache, falling back to small mapping (immune genes only)")
        ensg_map_path = os.path.join(os.path.dirname(os.path.dirname(BASE)),
                                     "output", "tcga_pancan", "ensg_map.json")
        with open(ensg_map_path) as f:
            sym2ensg = json.load(f)
        ensg2sym = {v: k for k, v in sym2ensg.items() if v}
    de["symbol"] = de["gene"].map(lambda g: ensg2sym.get(g, g))
    n_mapped = (de["symbol"] != de["gene"]).sum()
    print(f"symbol mapping: {n_mapped}/{len(de)} genes mapped successfully ({n_mapped/len(de)*100:.1f}%)")

    # 5. prerank GSEA
    rnk = de[["symbol", "z"]].drop_duplicates("symbol").sort_values("z", ascending=False)
    rnk = rnk.dropna()
    print(f"Ranked list: {len(rnk)} genes")
    rnk.to_csv(os.path.join(OUT, "fl_vs_ri_ranklist.csv"), index=False)

    gene_sets = [HALLMARK]
    try:
        gp.get_library(C7_IMMUNE, organism="Human")
        gene_sets.append(C7_IMMUNE)
    except Exception:
        print(f"  {C7_IMMUNE} unavailable, using only Hallmark")

    print("Running prerank GSEA ...")
    res = gp.prerank(rnk=rnk, gene_sets=gene_sets,
                     outdir=os.path.join(OUT, "gsea_out"),
                     min_size=5, max_size=500, permutation_num=1000,
                     seed=42, threads=4, no_plot=True)
    # Summary results (gseapy 1.3.1 output filename: gseapy.gene_set.prerank.report.csv)
    report = os.path.join(OUT, "gsea_out", "gseapy.gene_set.prerank.report.csv")
    if os.path.exists(report):
        t = pd.read_csv(report)
        t.to_csv(os.path.join(OUT, "gsea_summary.csv"), index=False)
        # Term column contains "library__pathway", extract pathway name
        t["Pathway"] = t["Term"].str.split("__", expand=True)[1].fillna(t["Term"])
        sig = t[t["FWER p-val"] < 0.25].sort_values("NES", ascending=False)
        print(f"\nEnrichment results: {len(t)} terms | FWER<0.05: {(t['FWER p-val'] < 0.05).sum()} terms | FWER<0.25: {len(sig)} terms")
        print("\n=== FL_high up-regulated pathways (NES>0) ===")
        up = sig[sig["NES"] > 0][["Pathway", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]]
        print(up.head(12).to_string(index=False))
        print("\n=== FL_high down-regulated pathways (NES<0, RI_high enriched) ===")
        dn = sig[sig["NES"] < 0][["Pathway", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]]
        print(dn.head(12).to_string(index=False))
    else:
        print(f"!! GSEA report not found: {report}")
        print("    Actual output:", os.listdir(os.path.join(OUT, "gsea_out")) if os.path.exists(os.path.join(OUT, "gsea_out")) else "no directory")

    print(f"\nOutput directory: {OUT}")

    # ---- Plot: significant pathway NES bar chart ----
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
            print(f"Figure saved: fig_gsea_fl_vs_ri.png ({len(sig)} significant pathways)")


if __name__ == "__main__":
    main()
