#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 Freeze ③ — GSEA headline pathways (freeze_a3_gsea.py)
===================================================
Purpose: extract from the gseapy.prerank raw output gsea_summary.csv the manuscript-cited
      6 headline pathways (FL-high positive: TNF-α/NF-κB, Inflammatory;
      FL-high negative / RI-high positive: E2F, G2-M, Myc, DNA repair),
      freeze into a3_gsea_headline.csv, for citation in Fig3 and the main text.

Input: article3/results/zp3_gsea_results/gsea_summary.csv (complete table of 50 pathways)
Output: article3/results/a3_gsea_headline.csv

Validation: checked against v0.2 manuscript
  TNF-α/NF-κB NES=+2.45 FDR<0.001；Inflammatory NES=+1.86；
  E2F NES=-2.12；G2-M NES=-2.10；Myc NES=-2.10；DNA repair NES=-1.67
"""
import os
import sys
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 levels = project root (verified empirically)
GSEA_CSV = os.path.join(ROOT, "article3", "results", "zp3_gsea_results", "gsea_summary.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_gsea_headline.csv")

# (keyword, manuscript NES, direction)
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
    print("=== GSEA headline pathway freeze ===")
    for kw, nes_exp, direction in HEADLINES:
        hit = g[g["Term"].str.contains(kw, case=False, regex=False)]
        if hit.empty:
            print(f"  !! {kw}: no hit"); ok_all = False; continue
        best = hit.sort_values("FDR q-val").iloc[0]
        nes, fdr = float(best["NES"]), float(best["FDR q-val"])
        ok_nes = abs(nes - nes_exp) < 0.05
        # Manuscript only declares FDR q<0.001 at TNF-α; actual for Inflammatory/DNA repair
        # FDR is 4.3e-3 / 4.9e-3 (q<0.01, not declared q<0.001 in manuscript) — assertion distinguishes.
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
        print(f"  {kw:<14} NES={nes:+.3f} (manuscript {nes_exp:+.2f}) FDR={fdr:.1e} "
              f"{'PASS' if ok else 'FAIL'}")
    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False)
    print(f"\nFrozen table: {OUT_CSV} ({len(res)} pathways)")
    print("RESULT:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
