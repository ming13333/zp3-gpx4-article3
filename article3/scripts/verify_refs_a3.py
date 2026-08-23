#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
verify_refs_a3.py — A3 参考文献核验 v2（Crossref DOI 解析 + PubMed PMID→DOI 反查交叉）
================================================================
方法（金标准双通道）：
  1. DOI 通道: Crossref works API 解析 → DOI 存在性 + 年份 + 期刊
  2. PMID 通道: NCBI eutils esummary 反查 PMID 的 DOI → 与稿件 DOI 比对
  状态:
    VERIFIED        : DOI 解析成功 + PMID 反查 DOI 与稿件一致
    VERIFIED_DOI    : DOI 解析成功（无 PMID 或 PMID 无 DOI 记录）
    SUSPICIOUS      : DOI 解析成功但 PMID→DOI 与稿件 DOI 不一致（需人工复核）
    FAIL            : DOI 无法解析（Crossref 无记录）
    NO_DOI          : 专著（74/75，无 DOI 正常）
输出: article3/results/ref_verify_report_a3.csv
"""
import os
import csv
import json
import time
import re
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))
MS = os.path.join(ROOT, "article3", "manuscripts", "Article3_A3英文初稿_v0.3.3.md")
OUT = os.path.join(ROOT, "article3", "results", "ref_verify_report_a3.csv")
MAILTO = "lm962272@gmail.com"
UA = f"a3-ref-verify/2.0 (mailto:{MAILTO})"


def extract_refs(text):
    m = re.search(r"## References\n(.*)$", text, re.S)
    refs = []
    for line in m.group(1).splitlines():
        line = line.strip()
        mm = re.match(r"^(\d+)\.\s+(.*)$", line)
        if mm:
            refs.append((int(mm.group(1)), mm.group(2)))
    return refs


def parse_line(line):
    doi = None
    md = re.search(r"doi:\s*(10\.\S+)", line, re.I)
    if md:
        doi = md.group(1).rstrip(".,;)")
    pmid = None
    mp = re.search(r"PMID:\s*(\d+)", line, re.I)
    if mp:
        pmid = mp.group(1)
    my = re.search(r"\b(19|20)\d{2}\b", line)
    year = my.group(0) if my else None
    # 本地标题：DOI/PMID 之前的完整文本（去作者段第一句）
    local = re.sub(r"\s*doi:\s*\S+\s*$", "", line, flags=re.I)
    local = re.sub(r"\s*PMID:\s*\d+\s*$", "", local, flags=re.I)
    return {"doi": doi, "pmid": pmid, "year": year, "local": local}


def cr_works(doi):
    url = f"https://api.crossref.org/works/{urllib.request.quote(doi)}?mailto={MAILTO}"
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            msg = json.loads(r.read().decode())["message"]
        title = (msg.get("title") or [""])[0]
        year = None
        for k in ("published-print", "published-online", "issued"):
            if msg.get(k) and msg[k].get("date-parts"):
                year = msg[k]["date-parts"][0][0]
                break
        journal = (msg.get("container-title") or [""])[0]
        return {"ok": True, "title": title, "year": year, "journal": journal,
                "vol": msg.get("volume", ""), "page": msg.get("page", "")}
    except Exception as e:
        return {"ok": False, "error": str(e)[:80]}


def pmid_doi(pmid):
    """NCBI eutils esummary 反查 PMID → DOI"""
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
           f"db=pubmed&id={pmid}&retmode=json&tool=refverify&email={MAILTO}")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read().decode())
        res = d["result"][pmid]
        for aid in res.get("articleids", []):
            if aid.get("idtype") == "doi":
                return aid.get("value")
        return None
    except Exception as e:
        return {"error": str(e)[:80]}


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", (s or "").lower()).strip()


def main():
    with open(MS, encoding="utf-8") as f:
        text = f.read()
    refs = extract_refs(text)
    print(f"提取参考文献: {len(refs)} 条")

    rows = []
    for num, line in refs:
        info = parse_line(line)
        base = {"num": num, "line": line[:80], "doi": info["doi"] or "",
                "pmid": info["pmid"] or "", "year_ms": info["year"] or ""}
        if not info["doi"]:
            rows.append({**base, "status": "NO_DOI", "cr_title": "", "cr_year": "",
                         "cr_journal": "", "pmid_doi": "", "note": "专著（无 DOI，正常）"})
            continue
        cr = cr_works(info["doi"])
        time.sleep(0.25)
        if not cr["ok"]:
            rows.append({**base, "status": "FAIL", "cr_title": "", "cr_year": "",
                         "cr_journal": "", "pmid_doi": "", "note": f"Crossref: {cr['error']}"})
            continue
        # PMID 反查
        pdoi = None
        pnote = ""
        if info["pmid"]:
            r = pmid_doi(info["pmid"])
            time.sleep(0.25)
            if isinstance(r, dict):
                pnote = f"eutils: {r['error']}"
            else:
                pdoi = r
        cr_doi_n = norm(info["doi"])
        pm_doi_n = norm(pdoi or "")
        if pdoi and pm_doi_n and pm_doi_n != cr_doi_n:
            status = "SUSPICIOUS"
            note = f"PMID→DOI({pdoi}) ≠ 稿件 DOI"
        elif pdoi:
            status = "VERIFIED"
            note = "DOI + PMID 双向一致"
        else:
            status = "VERIFIED_DOI"
            note = "DOI 解析成功（无 PMID 或 PMID 无 DOI）" + (f"; {pnote}" if pnote else "")
        # 年份对照（±1 容忍）
        yr_note = ""
        if info["year"] and cr["year"]:
            if abs(int(info["year"]) - int(cr["year"])) > 1:
                yr_note = f"; 年份: 稿件 {info['year']} vs Crossref {cr['year']}"
        rows.append({**base, "status": status, "cr_title": cr["title"][:90],
                     "cr_year": str(cr["year"] or ""), "cr_journal": cr["journal"][:50],
                     "pmid_doi": pdoi or "", "note": note + yr_note})

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    from collections import Counter
    cnt = Counter(r["status"] for r in rows)
    print(f"\n=== 结果 ===")
    for k in ("VERIFIED", "VERIFIED_DOI", "SUSPICIOUS", "FAIL", "NO_DOI"):
        print(f"  {k:14s}: {cnt.get(k, 0)}")
    print("\n--- SUSPICIOUS / FAIL 明细（需人工复核）---")
    for r in rows:
        if r["status"] in ("SUSPICIOUS", "FAIL"):
            print(f"  [{r['status']}] #{r['num']} {r['doi']} | {r['note']}")
    print(f"\n冻结: {OUT}")
    return 0 if cnt.get("FAIL", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())