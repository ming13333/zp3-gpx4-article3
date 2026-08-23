#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
probe_cptac_isoform_availability — 复现性探测：CPTAC-3 是否存在公开 mRNA isoform 级数据
================================================================
实测日期：2026-08-18（结论已写入 A3_isoform外部验证路线图.md）

背景：GDC CPTAC 页面 Data Types 表显示 "Isoform Expression Quantification TSV Open"，
需确认它是否为 mRNA 转录本（ENST）isoform，还是 miRNA isomiR。

方法：GDC Files API 两查询
  Q1: CPTAC-3 全部文件按 data_type 分面统计（facets=data_type, size=0）
  Q2: data_type=Isoform Expression Quantification 的文件名样本（验证 .mirnaseq. 前缀）
判定：
  - 若 Q2 文件名含 ".mirnaseq.isoforms." → 该 open 数据类型为 miRNA isomiR，非 mRNA isoform
  - 结合 Q1：
      * Gene Expression Quantification (open, 基因级)
      * Splice Junction Quantification (controlled, 需 dbGaP)
      * Isoform EQ == miRNA EQ 文件数 → isomiR 确认
  结论：CPTAC-3 无公开 mRNA 转录本 isoform 级数据

历史注记：2026-08-18 早前报告「GDC 无 isoform 级公开数据」部分正确（针对 mRNA）；
本脚本给出可复现证据链。原始数据源：https://gdc.cancer.gov/node/1179
"""
import json
import urllib.request
import urllib.parse

BASE = "https://api.gdc.cancer.gov/files"


def gdc_get(params):
    url = BASE + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=90) as r:
        return json.loads(r.read().decode())


def main():
    # Q1: CPTAC-3 data_type 分面
    f1 = {"op": "in", "content": {"field": "cases.project.project_id",
                                  "value": ["CPTAC-3"]}}
    q1 = gdc_get({"filters": json.dumps(f1), "facets": "data_type",
                  "format": "json", "size": 0})
    buckets = q1["data"]["aggregations"]["data_type"]["buckets"]
    out = {}
    for b in buckets:
        out[b["key"]] = b["doc_count"]
    print("=== CPTAC-3 全部文件 data_type 分布 ===")
    for k in sorted(out, key=lambda x: -out[x]):
        print(f"  {k:42s} {out[k]}")

    # Q2: Isoform Expression Quantification 文件样例（验证 mirnaseq 前缀）
    f2 = {"op": "and",
          "content": [
              {"op": "in", "content": {"field": "cases.project.project_id",
                                       "value": ["CPTAC-3"]}},
              {"op": "in", "content": {"field": "data_type",
                                       "value": ["Isoform Expression Quantification"]}},
          ]}
    q2 = gdc_get({"filters": json.dumps(f2), "fields": "file_name,access",
                  "format": "json", "size": 5})
    hits = q2["data"]["hits"]
    print("\n=== Isoform EQ 文件样例 ===")
    for h in hits:
        print(f"  [{h['access']}] {h['file_name']}")

    # 判定
    mirna_iso = sum(1 for h in hits if "mirnaseq" in h["file_name"])
    n_iso = out.get("Isoform Expression Quantification", 0)
    n_mir = out.get("miRNA Expression Quantification", 0)
    print(f"\nIsoform EQ 文件数 = {n_iso} | miRNA EQ 文件数 = {n_mir} | 相等? {n_iso == n_mir}")
    print(f"样例中 mirnaseq 前缀比例 = {mirna_iso}/{len(hits)}")
    if n_iso == n_mir and mirna_iso >= len(hits):
        print("结论: CPTAC-3 的 open 'Isoform Expression Quantification' = miRNA isomiR "
              "（非 mRNA 转录本 isoform）。mRNA isoform 级外部数据仍需 dbGaP 受控访问。")
        ok = True
    else:
        print("结论: 需人工复核，当前证据不足")
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())