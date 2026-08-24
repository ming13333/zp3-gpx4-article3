#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
probe_cptac_isoform_availability — Reproducibility probe: whether CPTAC-3 has public mRNA isoform-level data
================================================================
Test date: 2026-08-18 (conclusion written into A3_isoform_external_validation_roadmap.md)

Background: GDC CPTAC page Data Types table shows "Isoform Expression Quantification TSV Open",
Need to confirm whether it is an mRNA transcript (ENST) isoform, or miRNA isomiR.

Method: two queries to GDC Files API
  Q1: CPTAC-3 all files faceted by data_type (facets=data_type, size=0)
  Q2: sample filenames for data_type=Isoform Expression Quantification (verify .mirnaseq. prefix)
Determination:
  - If Q2 filenames contain ".mirnaseq.isoforms." → this open data type is miRNA isomiR, not mRNA isoform
  - Combined with Q1:
      * Gene Expression Quantification (open, gene-level)
      * Splice Junction Quantification (controlled, requires dbGaP)
      * Isoform EQ == miRNA EQ file count → isomiR confirmed
  Conclusion: CPTAC-3 has no public mRNA transcript isoform-level data

Historical note: On 2026-08-18, an earlier report that 'GDC has no isoform-level public data' was partially correct (for mRNA);
This script provides a reproducible evidence chain. Original data source: https://gdc.cancer.gov/node/1179
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
    # Q1: CPTAC-3 data_type faceting
    f1 = {"op": "in", "content": {"field": "cases.project.project_id",
                                  "value": ["CPTAC-3"]}}
    q1 = gdc_get({"filters": json.dumps(f1), "facets": "data_type",
                  "format": "json", "size": 0})
    buckets = q1["data"]["aggregations"]["data_type"]["buckets"]
    out = {}
    for b in buckets:
        out[b["key"]] = b["doc_count"]
    print("=== CPTAC-3 all files data_type distribution ===")
    for k in sorted(out, key=lambda x: -out[x]):
        print(f"  {k:42s} {out[k]}")

    # Q2: Isoform Expression Quantification file sample (verify mirnaseq prefix)
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
    print("\n=== Isoform EQ file sample ===")
    for h in hits:
        print(f"  [{h['access']}] {h['file_name']}")

    # Determination
    mirna_iso = sum(1 for h in hits if "mirnaseq" in h["file_name"])
    n_iso = out.get("Isoform Expression Quantification", 0)
    n_mir = out.get("miRNA Expression Quantification", 0)
    print(f"\nIsoform EQ file count = {n_iso} | miRNA EQ file count = {n_mir} | equal? {n_iso == n_mir}")
    print(f"sample mirnaseq prefix ratio = {mirna_iso}/{len(hits)}")
    if n_iso == n_mir and mirna_iso >= len(hits):
        print("Conclusion: CPTAC-3 open 'Isoform Expression Quantification' = miRNA isomiR "
              "(not mRNA transcript isoforms). mRNA isoform-level external data still requires dbGaP controlled access.")
        ok = True
    else:
        print("Conclusion: manual review needed, current evidence insufficient")
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
