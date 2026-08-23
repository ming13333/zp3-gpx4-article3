#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结④ — Mixed-effects isoform-immune 关联（freeze_a3_mixed_model.py）
========================================================================
目的：从长表（psi_immune_joined_samples.csv，9186 样本 × 32 癌种）独立复算
      mixed-effects 模型（未调整 + 调整 ZP3 总表达），冻结为
      a3_mixed_model_frozen.csv，供 Fig4 与正文引用。

模型（与 zp3_mixed_model.py / zp3_mixed_model_adjusted.py 一致）：
  - 未调整: score ~ PSI + (1 | Cancer)
  - 调整后: score ~ PSI + ZP3_total + (1 | Cancer)
  - 两转录本（FL / RI）× 7 免疫特征 = 14 行 × 2 组
  - BH FDR 校正（手工实现，与原脚本一致）

输入：
  - article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv
  - output/phase1_knowledge_gap_filling/gpx4_zp3_expr_matrix.csv（ZP3 总表达，按样本对齐）
输出：
  - article3/results/a3_mixed_model_frozen.csv
"""
import os
import sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根（实测验证）
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
    print(f"合并后: {len(df)} 样本 × {df['Cancer'].nunique()} 癌种（≥30 样本过滤后）")

    unadj = fit_models(df, "unadjusted")
    adj = fit_models(df, "adjusted_ZP3_total", add_total=True)
    out = pd.concat([unadj, adj], ignore_index=True)
    out.to_csv(OUT_CSV, index=False)

    print(f"\n冻结表写入: {OUT_CSV} ({len(out)} 行)")
    print("\n=== 关键数值（Fig4/正文引用）===")
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
    print(f"\n未调整: FL 正显著 {n_fl}/7 | RI 负显著 {n_ri}/7")
    print(f"调整总表达后: FL 正显著 {n_fla}/7")

    # 与稿件核对
    checks = {
        "M2 unadj": (fl[fl["Feature"] == "M2_Macrophage"], 0.2822, 8.44e-35),
        "M2 adj": (fla[fla["Feature"] == "M2_Macrophage"], 0.2383, 8.14e-24),
        "Cytolytic unadj.p": (fl[fl["Feature"] == "Cytolytic_activity"], np.nan, 0.578),
    }
    ok_all = True
    print("\n=== 与 v0.2 稿件核对 ===")
    for label, (sub, beta_exp, p_exp) in checks.items():
        r = sub.iloc[0]
        ok_b = pd.isna(beta_exp) or abs(r["Coef"] - beta_exp) < 0.001
        ok_p = abs(np.log10(r["P"]) - np.log10(p_exp)) < 0.3
        ok = bool(ok_b) and bool(ok_p)
        ok_all &= ok
        print(f"  {label}: β={r['Coef']:.4f} (稿 {beta_exp}) p={r['P']:.2e} (稿 {p_exp:.1e}) "
              f"{'PASS' if ok else 'FAIL'}")
    if n_fl != 6 or n_ri != 6 or n_fla != 6:
        ok_all = False
        print(f"  !! 显著计数不符: FL {n_fl}/6, RI {n_ri}/6, adj FL {n_fla}/6")
    else:
        print(f"  FL/RI/adj-FL 均为 6/7 正/负显著 PASS")
    print("\nRESULT:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()