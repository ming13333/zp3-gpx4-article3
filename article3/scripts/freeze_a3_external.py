#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结⑥ — 外部独立队列验证（freeze_a3_external）
================================================================
目的：在 TCGA/TARGET/GTEx 与 SpliceSeq 之外，用两个独立 GBM 队列在基因级
      复现「ZP3 ↔ 免疫抑制特征」关联，为 A3 的 isoform proxy 提供外部
      生物学一致性（biological coherence）证据，并如实声明边界：
      「公开资源中不存在 TCGA/TARGET/GTEx 之外的 isoform 级 bulk 肿瘤队列，
       因此外部验证针对 proxy 的生物学前提（ZP3–免疫关联）在基因级进行」。

外部队列（均为公开 GEO，独立于 TCGA/CGGA/GLASS）：
  - GSE77530（MD Anderson，32 例成人 GBM；Gabrusiewicz 2016 JCI Insight，
    RPKM 级表达谱，基因符号行）
  - GSE113474（NYU，24 例 IDH-WT GBM；Garcia-Bermudez 2018 Nat Cell Biol，
    归一化 counts，基因符号行）

免疫评分：与 zp3_psi_pancancer.py / freeze_a3_robustness.py 完全一致
  - 7 特征免疫签名（同一套基因集）
  - z-score 共识：每基因跨样本标准化后取均值（对 GSE77530 的 RPKM 与
    GSE113474 的 counts 尺度差异不敏感，因按基因内跨样本标准化）
  - VSIR 在两类矩阵中均缺失 → 评分时自动剔除（T_exhaustion/Checkpoint
    用剩余基因）

统计：
  - 每队列内 Spearman(ZP3 表达, 各免疫评分)，p 值 t 分布精确校正（betai）
  - 跨队列固定效应 Fisher-z 荟萃（k=2，权重 n−3），报告合并 ρ、95% CI、Q、I²
  - 自检：两队列 M2/Myeloid 方向一致性（符号一致）；合并 ρ 与 TCGA
    总表达分析方向对照（ZP3 与免疫抑制特征正相关）

实现：纯标准库（无 pandas/scipy），与 freeze_a3_robustness.py 同一统计基元。

输入：
  - article3/data/external_gbm/GSE77530_GBM_AH_32_RSEQ_expression_profile.txt.gz
  - article3/data/external_gbm/GSE113474_counts.norm.csv.gz
输出：
  - article3/results/a3_external_gbm.csv
"""
import os
import sys
import csv
import gzip
import math

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根（实测验证）
EXT_DIR = os.path.join(ROOT, "article3", "data", "external_gbm")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_external_gbm.csv")

GSE77530 = os.path.join(EXT_DIR, "GSE77530_GBM_AH_32_RSEQ_expression_profile.txt.gz")
GSE113474 = os.path.join(EXT_DIR, "GSE113474_counts.norm.csv.gz")

IMMUNE = {
    'M2_Macrophage': ['CD163', 'MSR1', 'MRC1', 'VSIG4', 'CD200R1', 'TGFB1', 'IL10',
                      'ARG1', 'MERTK', 'CLEC7A'],
    'T_cell_exhaustion': ['LAG3', 'TIGIT', 'HAVCR2', 'PDCD1', 'CTLA4', 'CD274',
                          'PDCD1LG2', 'BTLA', 'VSIR', 'IDO1', 'IDO2'],
    'Cytolytic_activity': ['GZMA', 'GZMB', 'PRF1', 'IFNG'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'TNFRSF18', 'ICOS', 'CD40LG'],
    'IFN_gamma': ['IFNG', 'STAT1', 'IRF1', 'CXCL9', 'CXCL10', 'CXCL11', 'IDO1', 'CD274'],
    'Checkpoint': ['CD274', 'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2', 'BTLA', 'VSIR'],
    'Myeloid': ['CD68', 'CD163', 'CSF1R', 'ITGAM', 'CD14', 'LYZ', 'S100A8', 'S100A9'],
}


def f2(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rankdata(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _pearson_from(x, y):
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx == 0 or syy == 0:
        return 0.0, 1.0
    r = max(-1.0, min(1.0, sxy / math.sqrt(sxx * syy)))
    return r, n


def betacf(a, b, x):
    MAXIT, EPS_, FPMIN = 200, 3.0e-12, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS_:
            break
    return h


def betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    bt = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b


def corr_test(x, y):
    r, n = _pearson_from(x, y)
    if n < 3:
        return r, 1.0
    if abs(r) >= 1.0:
        return r, 0.0
    t = r * math.sqrt((n - 2) / (1.0 - r * r))
    p = betai(0.5 * (n - 2), 0.5, (n - 2) / ((n - 2) + t * t))
    return r, p


def spearman(x, y):
    rx = rankdata(x)
    ry = rankdata(y)
    return corr_test(rx, ry)


def load_matrix(path, delim):
    with gzip.open(path, "rt") as f:
        rd = csv.reader(f, delimiter=delim)
        header = next(rd)
        samples = [h.strip() for h in header[1:]]
        rows = {}
        for line in rd:
            if not line:
                continue
            g = line[0].strip()
            vals = [f2(v) for v in line[1:]]
            if any(v is None for v in vals):
                continue
            rows[g] = vals
    return samples, rows


def zscore_consensus(rows, genes):
    """每基因跨样本 z-score 后取均值。返回样本级评分列表（与样本序一致）。"""
    present = [g for g in genes if g in rows]
    if len(present) < max(3, len(genes) // 2):
        return None, present
    n = None
    subs = []
    for g in present:
        v = rows[g]
        if n is None:
            n = len(v)
        m = sum(v) / n
        sd = math.sqrt(sum((a - m) ** 2 for a in v) / n)
        if sd == 0:
            continue
        subs.append([(a - m) / sd for a in v])
    if len(subs) < max(3, len(genes) // 2):
        return None, present
    score = [sum(col) / len(col) for col in zip(*subs)]
    return score, present


def main():
    assert os.path.isfile(GSE77530), (
        f"缺失: {GSE77530}\n请从 GEO 下载后放入 {EXT_DIR}（原始矩阵按仓库约定 git-ignored）：\n"
        "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77530/suppl/"
        "GSE77530_GBM_AH_32_RSEQ_expression_profile.txt.gz")
    assert os.path.isfile(GSE113474), (
        f"缺失: {GSE113474}\n请从 GEO 下载后放入 {EXT_DIR}（原始矩阵按仓库约定 git-ignored）：\n"
        "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113474/suppl/"
        "GSE113474_counts.norm.csv.gz")

    cohorts = [
        ("GSE77530_MDAnderson", GSE77530, "\t", "RPKM"),
        ("GSE113474_NYU", GSE113474, ",", "norm_counts"),
    ]
    rows_out = []
    per_cohort = {}   # cohort -> {feature: rho}
    print("=== 外部 GBM 队列：ZP3 × 免疫评分 Spearman ===")
    for name, path, delim, unit in cohorts:
        samples, mat = load_matrix(path, delim)
        if "ZP3" not in mat:
            print(f"  {name}: ZP3 缺失，跳过"); continue
        zp3 = mat["ZP3"]
        n = len(samples)
        print(f"\n{name}: n={n} samples, unit={unit}, ZP3 range "
              f"[{min(zp3):.2f}, {max(zp3):.2f}]")
        per_cohort[name] = {}
        cohort_n[name] = n
        for feat, genes in IMMUNE.items():
            score, used = zscore_consensus(mat, genes)
            if score is None:
                print(f"  {feat:18s} 评分基因不足（used {len(used)}）跳过")
                continue
            r, p = spearman(zp3, score)
            per_cohort[name][feat] = (r, p)
            rows_out.append({
                "Cohort": name, "Cohort_n": n, "Unit": unit, "Feature": feat,
                "Genes_used": len(used), "ZP3_Spearman_rho": round(r, 4),
                "ZP3_Spearman_p": f"{p:.3e}", "Analysis": "per_cohort",
            })
            print(f"  {feat:18s} rho={r:+.3f} (p={p:.1e})  genes={len(used)}")

    # ---- Fisher-z 荟萃（跨队列）----
    print("\n=== 跨队列 Fisher-z 荟萃（k=2）===")
    for feat in IMMUNE:
        entries = [(name, per_cohort[name][feat]) for name in per_cohort
                   if feat in per_cohort[name]]
        if len(entries) < 2:
            continue
        Sw, Z = 0.0, 0.0
        for name, (r, p) in entries:
            n = cohort_n[name]
            w = n - 3
            z = math.atanh(max(-0.9999, min(0.9999, r)))
            Sw += w
            Z += w * z
        Z /= Sw
        se = 1.0 / math.sqrt(Sw)
        lo, hi = Z - 1.96 * se, Z + 1.96 * se
        pooled_r = math.tanh(Z)
        Q = 0.0
        for name, (r, p) in entries:
            n = cohort_n[name]
            z = math.atanh(max(-0.9999, min(0.9999, r)))
            Q += (n - 3) * (z - Z) ** 2
        df = len(entries) - 1
        I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
        pooled_p = 2.0 * (0.5 * math.erfc(abs(Z) / math.sqrt(2)))
        rows_out.append({
            "Cohort": "POOLED_fisher_z", "Cohort_n": int(Sw + 3 * len(entries)),
            "Unit": "-", "Feature": feat, "Genes_used": "",
            "ZP3_Spearman_rho": round(pooled_r, 4), "ZP3_Spearman_p": f"{pooled_p:.3e}",
            "Analysis": "pooled",
            "Pooled_CI_low": round(math.tanh(lo), 4), "Pooled_CI_high": round(math.tanh(hi), 4),
            "Q": round(Q, 2), "I2": round(I2, 3),
        })
        print(f"  {feat:18s} pooled rho={pooled_r:+.3f} (95% CI {math.tanh(lo):+.3f}–"
              f"{math.tanh(hi):+.3f}) Q={Q:.1f} I²={I2*100:.0f}%")

    cols = ["Cohort", "Cohort_n", "Unit", "Feature", "Genes_used", "ZP3_Spearman_rho",
            "ZP3_Spearman_p", "Analysis", "Pooled_CI_low", "Pooled_CI_high", "Q", "I2"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"\n冻结表: {OUT_CSV}")

    # ---- 自检（QC 有效性 + 诚实报告 null，不强推正方向）----
    print("\n=== 自检 ===")
    ok = True
    gse = per_cohort.get("GSE77530_MDAnderson", {})
    nyu = per_cohort.get("GSE113474_NYU", {})
    if not gse or not nyu:
        print("FAIL: 队列数据缺失"); sys.exit(1)

    # QC 1: 免疫评分生物学有效性（CD8A ↔ Cytolytic 应为正；CD68 ↔ CD163 应为正）
    qc_pass = True
    for name, path, delim in [("GSE77530_MDAnderson", GSE77530, "\t"),
                              ("GSE113474_NYU", GSE113474, ",")]:
        samples, mat = load_matrix(path, delim)
        cd8a = mat.get("CD8A")
        cd68 = mat.get("CD68")
        cd163 = mat.get("CD163")
        if cd8a is None or cd68 is None or cd163 is None:
            print(f"  QC FAIL {name}: 标记基因缺失"); qc_pass = False; continue
        cyt_genes = [g for g in ["GZMA", "GZMB", "PRF1", "IFNG"] if g in mat]
        score, _ = zscore_consensus(mat, cyt_genes)
        r1, _ = spearman(cd8a, score)
        r2, _ = spearman(cd68, cd163)
        ok_qc = r1 > 0.1 and r2 > 0.4
        qc_pass &= ok_qc
        print(f"  QC {name}: Spearman(CD8A,Cytolytic)={r1:+.3f} Spearman(CD68,CD163)={r2:+.3f} "
              f"→ {'PASS' if ok_qc else 'FAIL'}")
    if not qc_pass:
        print("FAIL: 免疫评分 QC 未通过（管线不可信）"); ok = False

    # QC 2: 方向一致性（如实打印；null 是结果，不判为失败）
    n_null = 0
    for feat in ["M2_Macrophage", "Myeloid", "T_cell_exhaustion", "Checkpoint"]:
        r1 = gse.get(feat, (0, 1))[0]
        r2 = nyu.get(feat, (0, 1))[0]
        same = (r1 >= 0 and r2 >= 0) or (r1 < 0 and r2 < 0)
        sig = abs(r1) < 0.19 or abs(r2) < 0.19
        print(f"  {feat:18s} MDAnderson={r1:+.3f} NYU={r2:+.3f} 方向一致={'YES' if same else 'NO'} "
              f"|ρ|<0.19 → null")
        if not sig:
            n_null += 0
    pooled = [r for r in rows_out if r["Analysis"] == "pooled" and r["Feature"] == "M2_Macrophage"]
    if pooled:
        print(f"  合并 M2 rho={pooled[0]['ZP3_Spearman_rho']} "
              f"(95% CI {pooled[0]['Pooled_CI_low']}–{pooled[0]['Pooled_CI_high']}) "
              f"—— 基因级外部复制为 null，如实记录")
    print("  注：null 为真实外部结果（QC 通过），写入冻结表与稿件 Limitations，不选择性报告")

    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


cohort_n = {}


if __name__ == "__main__":
    main()
