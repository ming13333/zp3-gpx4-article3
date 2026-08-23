#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结③ — GSEA headline 通路（freeze_a3_gsea.py）
===================================================
目的：从 gseapy.prerank 原始输出 gsea_summary.csv 提取稿件引用的
      6 个 headline 通路（FL-high 正向：TNF-α/NF-κB、Inflammatory；
      FL-high 负向 / RI-high 正向：E2F、G2-M、Myc、DNA repair），
      冻结为 a3_gsea_headline.csv，供 Fig3 与正文引用。

输入：article3/results/zp3_gsea_results/gsea_summary.csv（50 通路整表）
输出：article3/results/a3_gsea_headline.csv

校验：与 v0.2 稿件核对
  TNF-α/NF-κB NES=+2.45 FDR<0.001；Inflammatory NES=+1.86；
  E2F NES=-2.12；G2-M NES=-2.10；Myc NES=-2.10；DNA repair NES=-1.67
"""
import os
import sys
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根（实测验证）
GSEA_CSV = os.path.join(ROOT, "article3", "results", "zp3_gsea_results", "gsea_summary.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_gsea_headline.csv")

# (关键词, 稿件 NES, 方向)
HEADLINES = [
    ("TNF-alpha", +2.45, "FL-high enriched"),
    ("Inflammatory", +1.86, "FL-high enriched"),
    ("E2F", -2.12, "RI-high enriched"),
    ("G2-M", -2.10, "RI-high enriched"),
    ("Myc", -2.10, "RI-high enriched"),
    ("DNA repair", -1.67, "RI-high enriched"),
]


def main():
    g = pd.read_csv(GSEA_CSV)
    rows = []
    ok_all = True
    print("=== GSEA headline 通路冻结 ===")
    for kw, nes_exp, direction in HEADLINES:
        hit = g[g["Term"].str.contains(kw, case=False, regex=False)]
        if hit.empty:
            print(f"  !! {kw}: 未命中"); ok_all = False; continue
        best = hit.sort_values("FDR q-val").iloc[0]
        nes, fdr = float(best["NES"]), float(best["FDR q-val"])
        ok_nes = abs(nes - nes_exp) < 0.05
        # 稿件仅在 TNF-α 处声明 FDR q<0.001；Inflammatory/DNA repair 的真实
        # FDR 为 4.3e-3 / 4.9e-3（q<0.01，未在稿件声明 q<0.001）——断言区分。
        strict_q = kw in ("TNF-alpha", "E2F", "G2-M", "Myc")
        ok_fdr = fdr < (0.001 if strict_q else 0.01)
        ok = ok_nes and ok_fdr
        ok_all &= ok
        rows.append({
            "Pathway": best["Term"], "Keyword": kw,
            "NES": round(nes, 4), "FDR_q": round(float(fdr), 6),
            "NOM_p": round(float(best["NOM p-val"]), 6),
            "ES": round(float(best["ES"]), 4),
            "Direction": direction, "Lead_genes": str(best["Lead_genes"])[:200],
        })
        print(f"  {kw:<14} NES={nes:+.3f} (稿 {nes_exp:+.2f}) FDR={fdr:.1e} "
              f"{'PASS' if ok else 'FAIL'}")
    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False)
    print(f"\n冻结表: {OUT_CSV} ({len(res)} 通路)")
    print("RESULT:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()