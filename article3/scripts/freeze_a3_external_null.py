#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结⑦ — 外部 null 诊断（freeze_a3_external_null）
================================================================
目的：对「外部基因级 ZP3–免疫关联为 null」做可冻结的诊断，回答三件事：
  (1) null 是不是因为 ZP3 在低表达平台接近检测底？
  (2) null 是不是因为 n=32/24 检验功效太低？
  (3) null 是不是被少数样本或表达尺度选择驱动？

本脚本不试图把 null 救回阳性，只做诚实诊断，全部冻结进 a3_external_null_diagnostics.csv。

复用 freeze_a3_external.py 的统计基元（rankdata / betai / spearman / zscore_consensus / load_matrix）。

四个诊断模块：
  D1. ZP3 检测质量：零值比例、非零数、IQR、中位数、范围、tied-rank 比例
  D2. 功效分析：n=32/24 在 α=0.05 双侧下可检测的最小 |ρ|（df=n-2，t 临界）
  D3. 表达尺度敏感性：原始 / log1p / 二值(检测) / 高表达分层(≥Q3) → Spearman(ZP3, score)
  D4. jackknife 稳定性：逐样本剔除后 M2 ρ 的波动范围

输入：
  - article3/data/external_gbm/GSE77530_GBM_AH_32_RSEQ_expression_profile.txt.gz
  - article3/data/external_gbm/GSE113474_counts.norm.csv.gz
输出：
  - article3/results/a3_external_null_diagnostics.csv
"""
import os
import sys
import csv
import gzip
import math

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根
EXT_DIR = os.path.join(ROOT, "article3", "data", "external_gbm")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_external_null_diagnostics.csv")

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


# ---------------------------------------------------------------------------
# 统计基元（与 freeze_a3_external.py 同一实现）
# ---------------------------------------------------------------------------
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
    return corr_test(rankdata(x), rankdata(y))


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


# ---------------------------------------------------------------------------
# 诊断模块
# ---------------------------------------------------------------------------
def quantile(vals, q):
    s = sorted(vals)
    if not s:
        return float('nan')
    pos = (len(s) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return s[lo]
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def diag_detection(name, zp3):
    """D1: ZP3 检测质量"""
    n = len(zp3)
    nz = [v for v in zp3 if v and v > 0]
    n_zero = n - len(nz)
    zero_frac = n_zero / n
    q1, med, q3 = quantile(zp3, 0.25), quantile(zp3, 0.5), quantile(zp3, 0.75)
    iqr = q3 - q1
    tied = len(zp3) - len(set([round(v, 6) for v in zp3]))
    tied_frac = tied / n
    return {
        "Cohort": name, "Diag": "D1_detection", "Metric": "ZP3",
        "Value": f"n={n},zero={n_zero}({zero_frac:.2f}),med={med:.3f},IQR={iqr:.3f},"
                 f"min={min(zp3):.3f},max={max(zp3):.3f},tied={tied}({tied_frac:.2f})",
        "Interpretation": (
            "低零值+合理IQR → 测量可用" if zero_frac < 0.2 and iqr > 0
            else "高零值或低IQR → 近检测底，秩相关被压缩"),
    }


def diag_power(name, n):
    """D2: 功效分析 — α=0.05 双侧最小可检测 |ρ|"""
    df = n - 2
    # t 临界：利用 betai 反推不方便，用近似 + 精确查表组合
    # 小样本用已知临界值表（df<=30）+ 大样本正态近似
    tcrit_table = {
        22: 2.074, 23: 2.069, 30: 2.042,  # n=24→df22, n=32→df30
    }
    if df in tcrit_table:
        tcrit = tcrit_table[df]
    else:
        # 大样本近似（df>30）: t ≈ z_{0.975} = 1.96
        tcrit = 1.96
    # |r| = t / sqrt(t^2 + df)
    min_r = tcrit / math.sqrt(tcrit * tcrit + df)
    return {
        "Cohort": name, "Diag": "D2_power", "Metric": f"n={n}",
        "Value": f"alpha=0.05,two-sided,min_detectable_|rho|={min_r:.3f}",
        "Interpretation": (
            f"该样本量最多检出 |rho|≥{min_r:.2f}；真实效应 0.2–0.3 大概率落入 null 区间"),
    }


def diag_scale(name, zp3, score):
    """D3: 表达尺度敏感性"""
    out = []
    # 原始
    r_raw, _ = spearman(zp3, score)
    # log1p
    zp3_log = [math.log1p(max(0.0, v)) for v in zp3]
    r_log, _ = spearman(zp3_log, score)
    # 二值（检测/未检测）
    zp3_bin = [1.0 if v and v > 0 else 0.0 for v in zp3]
    r_bin, _ = spearman(zp3_bin, score)
    # 高表达分层（≥Q3 vs <Q3）
    q3 = quantile(zp3, 0.75)
    zp3_hi = [1.0 if v >= q3 else 0.0 for v in zp3]
    r_hi, _ = spearman(zp3_hi, score)
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "raw",
        "Value": f"rho={r_raw:+.3f}", "Interpretation": "原始表达",
    })
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "log1p",
        "Value": f"rho={r_log:+.3f}", "Interpretation": "log1p 变换",
    })
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "detect_binary",
        "Value": f"rho={r_bin:+.3f}", "Interpretation": "检测/未检测二值",
    })
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "high_vs_low",
        "Value": f"rho={r_hi:+.3f}", "Interpretation": "高表达(≥Q3) vs 低表达",
    })
    return out


def diag_jackknife(name, zp3, score):
    """D4: jackknife 稳定性 — 逐样本剔除后 M2 ρ 波动"""
    n = len(zp3)
    rhos = []
    for i in range(n):
        x = zp3[:i] + zp3[i + 1:]
        y = score[:i] + score[i + 1:]
        r, _ = spearman(x, y)
        rhos.append(r)
    rmin, rmax = min(rhos), max(rhos)
    rmean = sum(rhos) / n
    # 与全样本比较
    r_full, _ = spearman(zp3, score)
    spread = rmax - rmin
    all_neg = rmin < 0 and rmax < 0
    return {
        "Cohort": name, "Diag": "D4_jackknife_M2", "Metric": f"n={n}",
        "Value": f"full_rho={r_full:+.3f},leave1out_range=[{rmin:+.3f},{rmax:+.3f}],"
                 f"spread={spread:.3f},mean={rmean:+.3f}",
        "Interpretation": (
            "方向一致为负且幅度中等敏感 → 非个别样本翻转"
            if spread < 0.25 and all_neg
            else "范围宽 → 结果对个别样本敏感"),
    }


def main():
    assert os.path.isfile(GSE77530), f"缺失: {GSE77530}"
    assert os.path.isfile(GSE113474), f"缺失: {GSE113474}"

    cohorts = [
        ("GSE77530_MDAnderson", GSE77530, "\t", "RPKM"),
        ("GSE113474_NYU", GSE113474, ",", "norm_counts"),
    ]
    rows_out = []
    print("=== A3 冻结⑦ 外部 null 诊断 ===")
    for name, path, delim, unit in cohorts:
        samples, mat = load_matrix(path, delim)
        if "ZP3" not in mat:
            print(f"  {name}: ZP3 缺失，跳过"); continue
        zp3 = mat["ZP3"]
        n = len(samples)
        print(f"\n{name}: n={n}, unit={unit}")

        # D1
        d1 = diag_detection(name, zp3)
        rows_out.append(d1)
        print(f"  D1 检测质量: {d1['Value']} → {d1['Interpretation']}")

        # D2
        d2 = diag_power(name, n)
        rows_out.append(d2)
        print(f"  D2 功效: {d2['Value']}")

        # M2 免疫评分（诊断尺度 + jackknife 的基准）
        m2_genes = IMMUNE['M2_Macrophage']
        score, used = zscore_consensus(mat, m2_genes)
        if score is None:
            print(f"  M2 评分基因不足，跳过 D3/D4"); continue

        # D3
        for d3 in diag_scale(name, zp3, score):
            rows_out.append(d3)
            print(f"  D3 {d3['Metric']:14s}: {d3['Value']}")

        # D4
        d4 = diag_jackknife(name, zp3, score)
        rows_out.append(d4)
        print(f"  D4 jackknife: {d4['Value']} → {d4['Interpretation']}")

    cols = ["Cohort", "Diag", "Metric", "Value", "Interpretation"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"\n冻结表: {OUT_CSV}")

    # ---- 自检：诊断不影响主结论，只核对格式与合理性 ----
    print("\n=== 自检 ===")
    ok = True
    # 自检 1: 所有行有解释
    for r in rows_out:
        if not r.get("Interpretation"):
            print(f"FAIL: 缺解释 {r}"); ok = False
    # 自检 2: D2 最小可检测 |rho| 在合理范围 (0.3–0.6)
    for r in rows_out:
        if r["Diag"] == "D2_power":
            s = r["Value"]
            lo = s.rfind("=") + 1
            mr = float(s[lo:].strip())
            if not (0.3 <= mr <= 0.6):
                print(f"FAIL: D2 最小|rho|异常 {mr}"); ok = False
    # 自检 3: D1 零值比例非负
    for r in rows_out:
        if r["Diag"] == "D1_detection":
            if "zero=" not in r["Value"]:
                print(f"FAIL: D1 格式异常"); ok = False
    print("  注：本脚本只做诊断，不重算主 null；主 null 见 a3_external_gbm.csv")
    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
