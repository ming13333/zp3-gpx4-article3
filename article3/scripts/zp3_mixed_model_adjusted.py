# -*- coding: utf-8 -*-
"""
A3 P1 补强: 调整 ZP3 总表达后的 isoform-immune 混合模型
=====================================================
核心主张检验: "isoform 信息超越 total expression"。
模型: score ~ FL_proportion + total_ZP3 + (1 | Cancer)
比较: 未调整模型(score ~ FL_proportion + (1|Cancer)) 与 调整后模型。
若 FL proportion 固定效应在控制总表达后仍显著 -> isoform 携带独立信息得到支持。

输入: psi_immune_joined_samples.csv + gpx4_zp3_expr_matrix.csv (按样本 join)
产物: zp3_psi_results/mixed_model_adjusted_results.csv
"""
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import rankdata

BASE = os.path.dirname(os.path.abspath(__file__))
JOINED = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_pancancer_results", "psi_immune_joined_samples.csv")
EXPR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "output", "phase1_knowledge_gap_filling", "gpx4_zp3_expr_matrix.csv")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_results")
os.makedirs(OUT, exist_ok=True)

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

def main():
    print("=== A3 P1: 调整总表达后的 mixed-effects ===\n")
    psi = pd.read_csv(JOINED)
    expr = pd.read_csv(EXPR)
    expr.columns = ["Sample"] + list(expr.columns[1:])
    df = psi.merge(expr[["Sample", "ZP3"]], on="Sample", how="inner")
    print(f"合并: {len(df)} 样本, {df['Cancer'].nunique()} 癌种")

    # 过滤癌种样本过少（<30），与原脚本一致
    sizes = df["Cancer"].value_counts()
    keep = sizes[sizes >= 30].index
    df = df[df["Cancer"].isin(keep)].copy()
    print(f"过滤后: {len(df)} 样本, {df['Cancer'].nunique()} 癌种\n")

    rows = []
    for tx, txlabel in [(FL, "FL canonical"), (RI, "Retained-intron")]:
        for feat in FEATURES:
            d = df[[tx, f"score_{feat}", "ZP3", "Cancer"]].dropna().copy()
            d.columns = ["psi", "score", "zp3", "Cancer"]
            if len(d) < 50:
                continue
            # 未调整模型（复现原始结论基线）
            try:
                m0 = smf.mixedlm("score ~ psi", d, groups=d["Cancer"]).fit()
                c0, p0 = m0.params["psi"], m0.pvalues["psi"]
                ci0 = m0.conf_int().loc["psi"]
            except Exception:
                c0, p0, ci0 = np.nan, np.nan, (np.nan, np.nan)
            # 调整总表达模型
            try:
                m1 = smf.mixedlm("score ~ psi + zp3", d, groups=d["Cancer"]).fit()
                c1, p1 = m1.params["psi"], m1.pvalues["psi"]
                ci1 = m1.conf_int().loc["psi"]
                c_zp3, p_zp3 = m1.params["zp3"], m1.pvalues["zp3"]
            except Exception as e:
                c1, p1, ci1 = np.nan, np.nan, (np.nan, np.nan)
                c_zp3, p_zp3 = np.nan, np.nan
                print(f"  ERR {tx} {feat}: {str(e)[:60]}")
            rows.append({
                "Transcript": tx, "Tx_Label": txlabel, "Feature": feat,
                "N": len(d), "N_cancer": d["Cancer"].nunique(),
                # unadjusted
                "Coef_unadj": round(float(c0), 4) if not np.isnan(c0) else np.nan,
                "P_unadj": float(p0) if not np.isnan(p0) else np.nan,
                "CI_unadj_low": round(float(ci0[0]), 4) if not np.isnan(c0) else np.nan,
                "CI_unadj_high": round(float(ci0[1]), 4) if not np.isnan(c0) else np.nan,
                # adjusted for total expression
                "Coef_adj": round(float(c1), 4) if not np.isnan(c1) else np.nan,
                "P_adj": float(p1) if not np.isnan(p1) else np.nan,
                "CI_adj_low": round(float(ci1[0]), 4) if not np.isnan(c1) else np.nan,
                "CI_adj_high": round(float(ci1[1]), 4) if not np.isnan(c1) else np.nan,
                "Coef_ZP3_total": round(float(c_zp3), 4) if not np.isnan(c_zp3) else np.nan,
                "P_ZP3_total": float(p_zp3) if not np.isnan(p_zp3) else np.nan,
            })

    res = pd.DataFrame(rows)
    # 对未调整和调整后分别算 FDR
    for col in ["P_unadj", "P_adj"]:
        sub = res[col].notna()
        if sub.sum() > 0:
            fdr = np.full(len(res), np.nan)
            fdr[sub.values] = bh_fdr(res.loc[sub, col].values)
            res[f"FDR_{col[2:]}"] = fdr  # FDR_unadj / FDR_adj

    res.to_csv(os.path.join(OUT, "mixed_model_adjusted_results.csv"), index=False)
    print(f"结果已存: mixed_model_adjusted_results.csv ({len(res)} 行)\n")

    print("=== FL canonical: 未调整 vs 调整总表达 ===")
    flr = res[res["Transcript"] == FL]
    for _, r in flr.iterrows():
        sig0 = "***" if r["FDR_unadj"] < 0.001 else ("**" if r["FDR_unadj"] < 0.01 else ("*" if r["FDR_unadj"] < 0.05 else ""))
        sig1 = "***" if r["FDR_adj"] < 0.001 else ("**" if r["FDR_adj"] < 0.01 else ("*" if r["FDR_adj"] < 0.05 else ""))
        print(f"  {r['Feature']:<22} unadj {r['Coef_unadj']:+.4f}(p={r['P_unadj']:.1e}{sig0}) "
              f"-> adj {r['Coef_adj']:+.4f}(p={r['P_adj']:.1e}{sig1}) | ZP3_total {r['Coef_ZP3_total']:+.4f}")
    n_fl_keep = ((flr["FDR_adj"] < 0.05) & (flr["Coef_adj"] > 0)).sum()
    print(f"\nFL 在调整总表达后仍正显著(FDR<0.05): {n_fl_keep}/7")

    print("\n=== Retained-intron: 未调整 vs 调整总表达 ===")
    rir = res[res["Transcript"] == RI]
    for _, r in rir.iterrows():
        sig0 = "***" if r["FDR_unadj"] < 0.001 else ("**" if r["FDR_unadj"] < 0.01 else ("*" if r["FDR_unadj"] < 0.05 else ""))
        sig1 = "***" if r["FDR_adj"] < 0.001 else ("**" if r["FDR_adj"] < 0.01 else ("*" if r["FDR_adj"] < 0.05 else ""))
        print(f"  {r['Feature']:<22} unadj {r['Coef_unadj']:+.4f}(p={r['P_unadj']:.1e}{sig0}) "
              f"-> adj {r['Coef_adj']:+.4f}(p={r['P_adj']:.1e}{sig1}) | ZP3_total {r['Coef_ZP3_total']:+.4f}")
    n_ri_keep = ((rir["FDR_adj"] < 0.05) & (rir["Coef_adj"] < 0)).sum()
    print(f"\nRI 在调整总表达后仍负显著(FDR<0.05): {n_ri_keep}/7")

if __name__ == "__main__":
    main()