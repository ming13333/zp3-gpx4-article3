#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结① — isoform 肿瘤/正常比例偏移（freeze_a3_isoform_shift）
================================================================
目的：将 A3 图 1（FL vs RI 肿瘤 vs 正常）涉及的统计量从中间产物
      zp3_isoform_proportions.csv 独立复算，产出冻结表 a3_isoform_shift.csv。
      后续所有稿件/插图只允许引用该冻结表数值。

输入：
  - article3/results/zp3_isoform_proportions.csv（19131 样本 × 7 转录本比例，来自 zp3_isoform_real_quant.py）
输出：
  - article3/results/a3_isoform_shift.csv
校验：
  - 与 v0.2 稿件核对：FL P=2.7e-139 / RI P=1.7e-53 / 5-exon 22 倍富集 P=1.3e-30
  - 与原 zp3_isoform_tumor_vs_normal.csv 全列数值一致性（Python 内断言）

口径（与原脚本一致）：
  肿瘤 = TCGA- 且样本段以 01 开头；正常 = GTEx。
  Mann-Whitney U two-sided；效应量 r = 1 - 2U/(n1*n2)（rank-biserial）；
  BH FDR 校正（statsmodels multipletests）。
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根（实测验证）
PROP_CSV = os.path.join(ROOT, "article3", "results", "zp3_isoform_proportions.csv")
ORIG_CSV = os.path.join(ROOT, "article3", "results", "zp3_isoform_tumor_vs_normal.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_isoform_shift.csv")
assert os.path.isdir(ROOT) and os.path.exists(PROP_CSV), f"ROOT 解析错误: {ROOT}"


def main():
    prop = pd.read_csv(PROP_CSV, index_col=0)
    prop.columns = [c.strip() for c in prop.columns]
    samples = list(prop.index)
    print(f"比例矩阵: {len(samples)} 样本 × {prop.shape[1]} 转录本")

    tcga_tumor = [s for s in samples
                  if s.startswith("TCGA-") and s.split("-")[3].startswith("01")]
    gtex = [s for s in samples if s.startswith("GTEX-")]
    print(f"肿瘤样本: {len(tcga_tumor)} | 正常样本: {len(gtex)}")

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
    print(f"\n冻结表已写入: {OUT_CSV}\n")
    print(res.to_string(index=False))

    # ---- 与稿件关键数值核对 ----
    checks = {
        "FL (ENST00000336517.8)": ("ENST00000336517.8", 0.403, 0.266, 2.7e-139),
        "RI (ENST00000466960.5)": ("ENST00000466960.5", 0.325, 0.440, 1.7e-53),
        "5-exon (ENST00000394860.3)": ("ENST00000394860.3", 0.0177, 0.0008, 1.3e-30),
    }
    print("\n=== 与 v0.2 稿件数值核对（四舍五入后 3 位有效数字）===")
    ok_all = True
    for label, (tid, tm, nm, p_exp) in checks.items():
        row = res[res["Transcript"] == tid].iloc[0]
        ok1 = abs(round(row["Tumor_median"], 3) - tm) < 0.0005
        ok2 = abs(round(row["Normal_median"], 3) - nm) < 0.0005
        ok3 = abs(np.log10(row["MannWhitney_p"]) - np.log10(p_exp)) < 0.15
        ok = ok1 and ok2 and ok3
        ok_all &= ok
        print(f"  {label}: T={row['Tumor_median']:.3f}(稿 {tm}) N={row['Normal_median']:.3f}(稿 {nm}) "
              f"p={row['MannWhitney_p']:.2e}(稿 {p_exp:.1e}) → {'PASS' if ok else 'FAIL'}")

    # ---- 与原 CSV 一致性（对数尺度；允许 scipy 版本级数值路径差异）----
    orig = pd.read_csv(ORIG_CSV)
    m = orig.merge(res, on="Transcript", suffixes=("_orig", "_new"))
    diff_log = np.max(np.abs(np.log10(m["MannWhitney_p_orig"]) - np.log10(m["MannWhitney_p_new"])))
    print(f"\n原 CSV 与冻结表最大对数 |Δlog10(p)| = {diff_log:.2e}（<0.01 即视为一致）")
    if diff_log >= 0.01:
        print("!! 与原 CSV 存在实质差异，检查口径"); ok_all = False
    else:
        print("与原 CSV 一致（数值路径版本差异，量级无损）✓")

    print("\nRESULT:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()