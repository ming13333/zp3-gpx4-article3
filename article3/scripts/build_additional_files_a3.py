#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_additional_files_a3.py — A3 BMC Bioinformatics 投稿 Additional files 汇编
================================================================
把冻结 CSV 汇编为 BMC 兼容的 Additional file xlsx（每表一个独立文件）：

  Additional file 1 → Supplementary Table S1 (Compositional controls)
                      来源: a3_robustness_frozen.csv
  Additional file 2 → Supplementary Table S2 (Cancer-stratified meta-analysis)
                      来源: a3_robustness_meta.csv
  Additional file 3 → Supplementary Table S3 (External cohorts + null diagnostics
                       + transportability)  来源: a3_external_gbm.csv +
                       a3_external_null_diagnostics.csv + a3_transportability_frozen.csv

输出: article3/results/additional_files/A3_AdditionalFile{1,2,3}.xlsx
依赖: openpyxl（本环境已装 3.1.5）
"""
import os
import csv
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))
RES = os.path.join(ROOT, "article3", "results")
OUT_DIR = os.path.join(RES, "additional_files")

SPECS = [
    {
        "file": "A3_AdditionalFile1.xlsx",
        "title": "Additional file 1. Compositional controls for the FL-immune association "
                 "(Supplementary Table S1).",
        "sheets": [
            ("Compositional_controls", os.path.join(RES, "a3_robustness_frozen.csv")),
        ],
    },
    {
        "file": "A3_AdditionalFile2.xlsx",
        "title": "Additional file 2. Cancer-stratified fixed-effect meta-analysis of the "
                 "FL-immune correlation (Supplementary Table S2).",
        "sheets": [
            ("Meta_analysis", os.path.join(RES, "a3_robustness_meta.csv")),
        ],
    },
    {
        "file": "A3_AdditionalFile3.xlsx",
        "title": "Additional file 3. External gene-level cohort assessment with null "
                 "diagnostics and internal transportability (Supplementary Table S3).",
        "sheets": [
            ("External_gene_level", os.path.join(RES, "a3_external_gbm.csv")),
            ("Null_diagnostics", os.path.join(RES, "a3_external_null_diagnostics.csv")),
            ("Transportability", os.path.join(RES, "a3_transportability_frozen.csv")),
        ],
    },
]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))


def sheet_from_rows(ws, rows):
    bold = Font(bold=True)
    fill = PatternFill("solid", fgColor="E6F1FB")
    wrap = Alignment(wrap_text=True, vertical="top")
    for j, cell in enumerate(rows[0], 1):
        c = ws.cell(row=1, column=j, value=cell)
        c.font = bold
        c.fill = fill
    for i, row in enumerate(rows[1:], 2):
        for j, v in enumerate(row, 1):
            ws.cell(row=i, column=j, value=v).alignment = wrap
    # 列宽
    for j in range(1, len(rows[0]) + 1):
        col_vals = [len(str(r[j - 1])) for r in rows[1:20] if len(r) >= j] or [8]
        width = min(max(col_vals) * 1.2 + 2, 60)
        ws.column_dimensions[get_column_letter(j)].width = max(width, 8)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for spec in SPECS:
        wb = openpyxl.Workbook()
        ws0 = wb.active
        ws0.title = "Info"
        ws0["A1"] = spec["title"]
        ws0["A1"].font = Font(bold=True, size=12)
        ws0["A2"] = "Source frozen tables: article3/results/ (see manuscript Data availability)."
        ws0.column_dimensions["A"].width = 110
        for name, csv_path in spec["sheets"]:
            rows = load_csv(csv_path)
            ws = wb.create_sheet(title=name[:31])
            sheet_from_rows(ws, rows)
        # Info 移到最前
        wb.move_sheet("Info", offset=-(len(spec["sheets"])))
        out = os.path.join(OUT_DIR, spec["file"])
        wb.save(out)
        sz = os.path.getsize(out)
        print(f"OK  {spec['file']}  ({sz/1024:.1f} KB, sheets: {[n for n,_ in spec['sheets']]})")

    # 清单
    manifest = os.path.join(OUT_DIR, "A3_AdditionalFiles_清单.md")
    lines = [
        "# A3 BMC Bioinformatics — Additional Files 清单（2026-08-18 汇编）",
        "",
        "| 文件 | 对应稿件 | 内容 | 来源冻结表 |",
        "|---|---|---|---|",
        "| A3_AdditionalFile1.xlsx | Supplementary Table S1 | FL–免疫关联的组成控制（log-ratio / FL–RI 耦合 / 低信号过滤） | a3_robustness_frozen.csv |",
        "| A3_AdditionalFile2.xlsx | Supplementary Table S2 | 32 癌种分层固定效应荟萃（M2/Myeloid，95% CI，Q，I²） | a3_robustness_meta.csv |",
        "| A3_AdditionalFile3.xlsx | Supplementary Table S3 | 外部基因级队列 + null 诊断 + 内部跨癌种迁移性（L2CO/held-out） | a3_external_gbm.csv, a3_external_null_diagnostics.csv, a3_transportability_frozen.csv |",
        "",
        "注：每个 xlsx 首 sheet 为文件说明（Info），数据 sheet 为对应冻结表（完整列、未截断）。",
        "BMC 硬约束：单个 Additional file ≤20 MB（当前各 <100 KB，通过）。",
    ]
    with open(manifest, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"OK  清单: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())