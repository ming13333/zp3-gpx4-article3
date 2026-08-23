# -*- coding: utf-8 -*-
"""
Article 3 补强③ — Mixed-effects 加固层内/层间结论
=================================================
目的：用随机截距混合模型控制癌种效应，检验 PSI 对免疫评分的固定效应是否仍显著，
从而区分"癌种间差异（生态学混杂）"与"癌种内 PSI 驱动"。

模型：immune_score ~ PSI + (1 | Cancer)   （statsmodels MixedLM）
  - FL_PSI / RI_PSI 分别建模
  - 若固定效应显著 → 层内驱动独立于癌种 → Simpson 悖论"层内驱动"结论加固
  - 同时报告固定效应系数（每单位 PSI 的免疫评分变化）

输入：zp3_psi_pancancer_results/psi_immune_joined_samples.csv（9186 样本长表）
产物：zp3_psi_results/mixed_model_results.csv + fig_mixed_model.png
"""
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
JOINED = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_pancancer_results", "psi_immune_joined_samples.csv")
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_results")
os.makedirs(OUT, exist_ok=True)

FL = "ENST00000336517.8"
RI = "ENST00000466960.5"
FEATURES = ["M2_Macrophage", "T_cell_exhaustion", "Cytolytic_activity",
            "Treg", "IFN_gamma", "Checkpoint", "Myeloid"]


def main():
    print("=== Article 3 补强③: Mixed-effects 模型 ===\n")
    df = pd.read_csv(JOINED)
    print(f"长表: {len(df)} 样本 × {df.shape[1]} 列 | 癌种数: {df['Cancer'].nunique()}")

    # 过滤癌种样本过少（<30）
    sizes = df["Cancer"].value_counts()
    keep = sizes[sizes >= 30].index
    df = df[df["Cancer"].isin(keep)].copy()
    print(f"过滤后: {len(df)} 样本, {df['Cancer'].nunique()} 癌种")
    print("癌种分布:", dict(df["Cancer"].value_counts().head(8)))

    rows = []
    for tx, txlabel in [(FL, "FL canonical"), (RI, "Retained-intron")]:
        for feat in FEATURES:
            d = df[[tx, f"score_{feat}", "Cancer"]].dropna().copy()
            d.columns = ["psi", "score", "Cancer"]
            if len(d) < 50:
                continue
            try:
                m = smf.mixedlm("score ~ psi", d, groups=d["Cancer"]).fit()
                coef = m.params["psi"]
                se = m.bse["psi"]
                pval = m.pvalues["psi"]
                ci = m.conf_int().loc["psi"]
                # 组内相关（ICC 近似：随机截距方差占比）
                vc = m.cov_re.iloc[0, 0]
                resid = m.scale
                icc = vc / (vc + resid) if (vc + resid) > 0 else np.nan
                rows.append({
                    "Transcript": tx, "Tx_Label": txlabel, "Feature": feat,
                    "N": len(d), "N_cancer": d["Cancer"].nunique(),
                    "Coef": round(float(coef), 4), "SE": round(float(se), 4),
                    "Z": round(float(coef / se), 3) if se != 0 else np.nan,
                    "P": float(pval), "CI_low": round(float(ci[0]), 4),
                    "CI_high": round(float(ci[1]), 4), "ICC": round(float(icc), 4),
                })
            except Exception as e:
                rows.append({"Transcript": tx, "Tx_Label": txlabel, "Feature": feat,
                             "N": len(d), "N_cancer": d["Cancer"].nunique(),
                             "Coef": np.nan, "SE": np.nan, "Z": np.nan,
                             "P": np.nan, "CI_low": np.nan, "CI_high": np.nan,
                             "ICC": np.nan, "Error": str(e)[:80]})

    res = pd.DataFrame(rows)
    # BH FDR
    pv = res["P"].values
    if len(pv):
        from scipy.stats import rankdata
        n = len(pv)
        order = np.argsort(pv)
        ranked = pv[order]
        q = ranked * n / (np.arange(1, n + 1))
        q = np.minimum.accumulate(q[::-1])[::-1]
        fdr = np.empty(n)
        fdr[order] = q
        res["FDR"] = np.minimum(fdr, 1.0)
    res.to_csv(os.path.join(OUT, "mixed_model_results.csv"), index=False)
    print(f"\n结果已存: mixed_model_results.csv ({len(res)} 行)")

    # 汇总展示
    print("\n=== FL canonical 固定效应 ===")
    flr = res[res["Transcript"] == FL]
    for _, r in flr.iterrows():
        sig = "***" if r["FDR"] < 0.001 else ("**" if r["FDR"] < 0.01 else ("*" if r["FDR"] < 0.05 else ""))
        print(f"  {r['Feature']:<22} coef={r['Coef']:+.4f} (95% CI {r['CI_low']:+.4f}~{r['CI_high']:+.4f}) "
              f"p={r['P']:.2e} FDR={r['FDR']:.2e} {sig}")
    print("\n=== Retained-intron 固定效应 ===")
    rir = res[res["Transcript"] == RI]
    for _, r in rir.iterrows():
        sig = "***" if r["FDR"] < 0.001 else ("**" if r["FDR"] < 0.01 else ("*" if r["FDR"] < 0.05 else ""))
        print(f"  {r['Feature']:<22} coef={r['Coef']:+.4f} (95% CI {r['CI_low']:+.4f}~{r['CI_high']:+.4f}) "
              f"p={r['P']:.2e} FDR={r['FDR']:.2e} {sig}")

    n_fl = ((res["Transcript"] == FL) & (res["FDR"] < 0.05) & (res["Coef"] > 0)).sum()
    n_ri = ((res["Transcript"] == RI) & (res["FDR"] < 0.05) & (res["Coef"] < 0)).sum()
    print(f"\n=== 判定 ===")
    print(f"FL 正显著(FDR<0.05): {n_fl}/7 | RI 负显著(FDR<0.05): {n_ri}/7")
    if n_fl >= 4:
        print("→ 层内驱动结论【加固】：控制癌种随机效应后，FL-PSI 固定效应仍显著正")
    elif n_fl >= 1:
        print("→ 层内驱动部分支持（部分特征显著）")
    else:
        print("→ 层内驱动未获支持：固定效应不显著，需 sensitivity analysis")

    # 图：森林图（固定效应系数）
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    for ax, tx, title in [(axes[0], FL, "FL canonical (9-exon)"),
                          (axes[1], RI, "Retained-intron")]:
        sub = res[res["Transcript"] == tx].dropna(subset=["Coef"])
        if len(sub) == 0:
            ax.set_title(f"{title}\n(no data)", fontsize=11)
            continue
        sub = sub.sort_values("Coef")
        y = np.arange(len(sub))
        colors = ["#C00000" if (r["FDR"] < 0.05 and r["Coef"] > 0)
                  else "#1D9E75" if (r["FDR"] < 0.05 and r["Coef"] < 0)
                  else "#999999" for _, r in sub.iterrows()]
        for i, (_, r) in enumerate(sub.iterrows()):
            ax.errorbar([r["Coef"]], [y[i]],
                        xerr=[[r["Coef"] - r["CI_low"]], [r["CI_high"] - r["Coef"]]],
                        fmt="o", ecolor=colors[i], color=colors[i], ms=7,
                        capsize=3, linewidth=1.5, zorder=3)
        ax.axvline(0, color="grey", lw=0.8, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels(sub["Feature"], fontsize=9)
        ax.set_xlabel("Fixed-effect coef (immune score per unit PSI)", fontsize=10)
        ax.set_title(f"{title}\nMixed model: score ~ PSI + (1|Cancer)", fontsize=11)
        ax.grid(axis="x", alpha=0.3)
    fig.suptitle("Article 3 补强③: 控制癌种随机效应后的 PSI 固定效应", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_mixed_model.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n图已存: fig_mixed_model.png")


if __name__ == "__main__":
    main()
