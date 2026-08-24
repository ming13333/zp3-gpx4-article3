#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 Freeze ④ — Mixed-effects isoform-immune association (freeze_a3_mixed_model.py)
========================================================================
Purpose: Independently recompute from the long table (psi_immune_joined_samples.csv, 9186 samples × 32 cancer types)
      mixed-effects model (unadjusted + adjusted for ZP3 total expression), frozen as
      a3_mixed_model_frozen.csv, for citation in Fig4 and the main text.

Model (consistent with zp3_mixed_model.py / zp3_mixed_model_adjusted.py):
  - Unadjusted: score ~ PSI + (1 | Cancer)
  - Adjusted: score ~ PSI + ZP3_total + (1 | Cancer)
  - Two transcripts (FL / RI) × 7 immune features = 14 rows × 2 groups
  - BH FDR correction (hand-implemented, consistent with original script)

Input:
  - article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv
  - output/phase1_knowledge_gap_filling/gpx4_zp3_expr_matrix.csv (ZP3 total expression, aligned by sample)
Output:
  - article3/results/a3_mixed_model_frozen.csv
"""
import os
import sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 levels = project root (verified)
JOINED = os.path.join(ROOT, "article3", "results",
                      "zp3_psi_pancancer_results", "psi_immune_joined_samples.csv")
EXPR = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling", "gpx4_zp3_expr_matrix.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_mixed_model_frozen.csv")

FL = "ENST00000336517.8"
RI = "ENST00000466960.5"
FEATURES = ["M2_Macrophage", "T_cell_exhaustion", "Cytolytic_activity",
            "Treg", "IFN_gamma", "Checkpoint", "Myeloid"]


def bh_fdr(pvals):
    pv = np.asarray(pvals, dtype=float)
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = q
    return np.minimum(fdr, 1.0)


def fit_models(df, label, add_total=False):
    rows = []
    for tx, txlabel in [(FL, "FL canonical"), (RI, "Retained-intron")]:
        for feat in FEATURES:
            cols = [tx, f"score_{feat}", "Cancer"] + (["ZP3"] if add_total else [])
            d = df[cols].dropna().copy()
            d.columns = ["psi", "score", "Cancer"] + (["zp3"] if add_total else [])
            if len(d) < 50:
                continue
            try:
                if add_total:
                    m = smf.mixedlm("score ~ psi + zp3", d, groups=d["Cancer"]).fit()
                else:
                    m = smf.mixedlm("score ~ psi", d, groups=d["Cancer"]).fit()
                coef = float(m.params["psi"])
                se = float(m.bse["psi"])
                pval = float(m.pvalues["psi"])
                ci = m.conf_int().loc["psi"].values
                vc = m.cov_re.iloc[0, 0]
                resid = m.scale
                icc = vc / (vc + resid) if (vc + resid) > 0 else np.nan
                rows.append({
                    "Model": label, "Transcript": tx, "Tx_Label": txlabel,
                    "Feature": feat, "N": len(d), "N_cancer": d["Cancer"].nunique(),
                    "Coef": round(coef, 4), "SE": round(se, 4),
                    "Z": round(coef / se, 3) if se != 0 else np.nan,
                    "P": pval, "CI_low": round(float(ci[0]), 4),
                    "CI_high": round(float(ci[1]), 4), "ICC": round(float(icc), 4),
                })
            except Exception as e:
                rows.append({"Model": label, "Transcript": tx, "Tx_Label": txlabel,
                             "Feature": feat, "N": len(d), "N_cancer": d["Cancer"].nunique(),
                             "Coef": np.nan, "SE": np.nan, "Z": np.nan, "P": np.nan,
                             "CI_low": np.nan, "CI_high": np.nan, "ICC": np.nan,
                             "Error": str(e)[:80]})
    res = pd.DataFrame(rows)
    fdr = bh_fdr(res["P"].fillna(1.0).values)
    res["FDR"] = fdr
    res.loc[res["P"].isna(), "FDR"] = np.nan
    return res


def main():
    psi = pd.read_csv(JOINED)
    expr = pd.read_csv(EXPR)
    expr.columns = ["Sample"] + list(expr.columns[1:])
    df = psi.merge(expr[["Sample", "ZP3"]], on="Sample", how="inner")
    sizes = df["Cancer"].value_counts()
    keep = sizes[sizes >= 30].index
    df = df[df["Cancer"].isin(keep)].copy()
    print(f"After merging: {len(df)} samples × {df['Cancer'].nunique()} cancer types (after >=30 sample filter)")

    unadj = fit_models(df, "unadjusted")
    adj = fit_models(df, "adjusted_ZP3_total", add_total=True)
    out = pd.concat([unadj, adj], ignore_index=True)
    out.to_csv(OUT_CSV, index=False)

    print(f"Frozen table written: {OUT_CSV} ({len(out)} rows)")
    print("\n=== Key values (Fig4/main text citation) ===")
    fl = out[(out["Transcript"] == FL) & (out["Model"] == "unadjusted")]
    fla = out[(out["Transcript"] == FL) & (out["Model"] == "adjusted_ZP3_total")]
    for _, r in fl.iterrows():
        sig = "***" if r["FDR"] < 0.001 else ("**" if r["FDR"] < 0.01 else "*")
        print(f"  FL {r['Feature']:<20} β={r['Coef']:+.4f} p={r['P']:.2e} FDR={r['FDR']:.1e} {sig}")
    n_fl = ((fl["FDR"] < 0.05) & (fl["Coef"] > 0)).sum()
    n_ria = 0
    ri = out[(out["Transcript"] == RI) & (out["Model"] == "unadjusted")]
    n_ri = ((ri["FDR"] < 0.05) & (ri["Coef"] < 0)).sum()
    n_fla = ((fla["FDR"] < 0.05) & (fla["Coef"] > 0)).sum()
    print(f"\nUnadjusted: FL positive significant {n_fl}/7 | RI negative significant {n_ri}/7")
    print(f"After adjusting for total expression: FL positive significant {n_fla}/7")

    # cross-check with manuscript
    checks = {
        "M2 unadj": (fl[fl["Feature"] == "M2_Macrophage"], 0.2822, 8.44e-35),
        "M2 adj": (fla[fla["Feature"] == "M2_Macrophage"], 0.2383, 8.14e-24),
        "Cytolytic unadj.p": (fl[fl["Feature"] == "Cytolytic_activity"], np.nan, 0.578),
    }
    ok_all = True
    print("\n=== Check against v0.2 manuscript ===")
    for label, (sub, beta_exp, p_exp) in checks.items():
        r = sub.iloc[0]
        ok_b = pd.isna(beta_exp) or abs(r["Coef"] - beta_exp) < 0.001
        ok_p = abs(np.log10(r["P"]) - np.log10(p_exp)) < 0.3
        ok = bool(ok_b) and bool(ok_p)
        ok_all &= ok
        print(f"  {label}: β={r['Coef']:.4f} (draft {beta_exp}) p={r['P']:.2e} (draft {p_exp:.1e}) "
              f"{'PASS' if ok else 'FAIL'}")
    if n_fl != 6 or n_ri != 6 or n_fla != 6:
        ok_all = False
        print(f"  !! Significant count mismatch: FL {n_fl}/6, RI {n_ri}/6, adj FL {n_fla}/6")
    else:
        print(f"  FL/RI/adj-FL are all 6/7 positive/negative significant PASS")
    print("\nRESULT:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
