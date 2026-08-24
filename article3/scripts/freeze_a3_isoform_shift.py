#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 Freeze① — isoform tumor/normal proportion shift (freeze_a3_isoform_shift)
================================================================
Purpose: Independently recompute the statistics involved in A3 Figure 1 (FL vs RI tumor vs normal) from the intermediate product
      zp3_isoform_proportions.csv, and produce the frozen table a3_isoform_shift.csv.
      All subsequent manuscripts/figures may only cite values from this frozen table.

Input:
  - article3/results/zp3_isoform_proportions.csv (19131 samples × 7 transcript proportions, from zp3_isoform_real_quant.py)
Output:
  - article3/results/a3_isoform_shift.csv
Validation:
  - Check against v0.2 manuscript: FL P=2.7e-139 / RI P=1.7e-53 / 5-exon 22-fold enrichment P=1.3e-30
  - Column-wise numerical consistency with the original zp3_isoform_tumor_vs_normal.csv (Python assertion)

Definitions (consistent with original script):
  Tumor = TCGA- and sample segment starts with 01; Normal = GTEx.
  Mann-Whitney U two-sided; effect size r = 1 - 2U/(n1*n2) (rank-biserial);
  BH FDR correction (statsmodels multipletests).
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 levels = project root (verified empirically)
PROP_CSV = os.path.join(ROOT, "article3", "results", "zp3_isoform_proportions.csv")
ORIG_CSV = os.path.join(ROOT, "article3", "results", "zp3_isoform_tumor_vs_normal.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_isoform_shift.csv")
assert os.path.isdir(ROOT) and os.path.exists(PROP_CSV), f"ROOT parsing error: {ROOT}"


def main():
    prop = pd.read_csv(PROP_CSV, index_col=0)
    prop.columns = [c.strip() for c in prop.columns]
    samples = list(prop.index)
    print(f"Proportion matrix: {len(samples)} samples × {prop.shape[1]} transcripts")

    tcga_tumor = [s for s in samples
                  if s.startswith("TCGA-") and s.split("-")[3].startswith("01")]
    gtex = [s for s in samples if s.startswith("GTEX-")]
    print(f"Tumor samples: {len(tcga_tumor)} | Normal samples: {len(gtex)}")

    rows = []
    for tid in prop.columns:
        t_vals = prop.loc[tcga_tumor, tid].dropna()
        n_vals = prop.loc[gtex, tid].dropna()
        if len(t_vals) < 10 or len(n_vals) < 10:
            continue
        u, p = stats.mannwhitneyu(t_vals, n_vals, alternative="two-sided")
        r = 1 - (2 * u) / (len(t_vals) * len(n_vals))
        rows.append({
            "Transcript": tid,
            "Tumor_median": float(t_vals.median()),
            "Normal_median": float(n_vals.median()),
            "Tumor_mean": float(t_vals.mean()),
            "Normal_mean": float(n_vals.mean()),
            "Tumor_Q1": float(t_vals.quantile(0.25)),
            "Tumor_Q3": float(t_vals.quantile(0.75)),
            "Normal_Q1": float(n_vals.quantile(0.25)),
            "Normal_Q3": float(n_vals.quantile(0.75)),
            "Tumor_n": len(t_vals),
            "Normal_n": len(n_vals),
            "MannWhitney_p": float(p),
            "Effect_r": float(r),
        })
    res = pd.DataFrame(rows)
    _, fdr, _, _ = multipletests(res["MannWhitney_p"], method="fdr_bh")
    res["FDR"] = fdr
    res["Tumor_over_Normal_ratio"] = res["Tumor_median"] / res["Normal_median"].replace(0, np.nan)
    res = res.sort_values("MannWhitney_p")
    res.to_csv(OUT_CSV, index=False)
    print(f"\nFrozen table written to: {OUT_CSV}\n")
    print(res.to_string(index=False))

    # ---- Cross-check against manuscript key values ----
    checks = {
        "FL (ENST00000336517.8)": ("ENST00000336517.8", 0.403, 0.266, 2.7e-139),
        "RI (ENST00000466960.5)": ("ENST00000466960.5", 0.325, 0.440, 1.7e-53),
        "5-exon (ENST00000394860.3)": ("ENST00000394860.3", 0.0177, 0.0008, 1.3e-30),
    }
    print("\n=== Cross-check with v0.2 manuscript values (rounded to 3 significant digits) ===")
    ok_all = True
    for label, (tid, tm, nm, p_exp) in checks.items():
        row = res[res["Transcript"] == tid].iloc[0]
        ok1 = abs(round(row["Tumor_median"], 3) - tm) < 0.0005
        ok2 = abs(round(row["Normal_median"], 3) - nm) < 0.0005
        ok3 = abs(np.log10(row["MannWhitney_p"]) - np.log10(p_exp)) < 0.15
        ok = ok1 and ok2 and ok3
        ok_all &= ok
        print(f"  {label}: T={row['Tumor_median']:.3f}(manuscript {tm}) N={row['Normal_median']:.3f}(manuscript {nm}) "
              f"p={row['MannWhitney_p']:.2e}(manuscript {p_exp:.1e}) → {'PASS' if ok else 'FAIL'}")

    # ---- Consistency with original CSV (log scale; allow scipy version-level numerical path differences) ----
    orig = pd.read_csv(ORIG_CSV)
    m = orig.merge(res, on="Transcript", suffixes=("_orig", "_new"))
    diff_log = np.max(np.abs(np.log10(m["MannWhitney_p_orig"]) - np.log10(m["MannWhitney_p_new"])))
    print(f"\nMax log difference between original CSV and frozen table |Δlog10(p)| = {diff_log:.2e} (<0.01 considered consistent)")
    if diff_log >= 0.01:
        print("!! Substantial difference from original CSV, check definitions"); ok_all = False
    else:
        print("Consistent with original CSV (numerical path version difference, no loss in magnitude) ✓")

    print("\nRESULT:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
