#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 Freeze ⑦ - External Null Diagnostics (freeze_a3_external_null)
================================================================
Purpose: Provide freezable diagnostics for the external gene-level ZP3–immune association being null, answering three things:
  (1) Is the null result because ZP3 approaches the detection floor on a low-expression platform?
  (2) Is the null result because n=32/24 has too low statistical power?
  (3) Is the null result driven by a few samples or expression-scale choices?

This script does not attempt to rescue the null into positive; it only performs honest diagnostics, all frozen into a3_external_null_diagnostics.csv.

Reuses statistical primitives from freeze_a3_external.py (rankdata / betai / spearman / zscore_consensus / load_matrix).

Four diagnostic modules:
  D1. ZP3 detection quality: zero proportion, nonzero count, IQR, median, range, tied-rank proportion
  D2. Power analysis: minimum detectable |ρ| for n=32/24 under two-sided α=0.05 (df=n-2, t critical)
  D3. Expression-scale sensitivity: raw / log1p / binary (detection) / high-expression stratification (≥Q3) → Spearman(ZP3, score)
  D4. Jackknife stability: fluctuation range of M2 ρ after removing each sample one by one

Input:
  - article3/data/external_gbm/GSE77530_GBM_AH_32_RSEQ_expression_profile.txt.gz
  - article3/data/external_gbm/GSE113474_counts.norm.csv.gz
Output:
  - article3/results/a3_external_null_diagnostics.csv
"""
import os
import sys
import csv
import gzip
import math

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 levels = project root
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
# Statistical primitives (same implementation as freeze_a3_external.py)
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
# Diagnostics module
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
    """D1: ZP3 detection quality"""
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
            "Low zero fraction + reasonable IQR → measurement usable" if zero_frac < 0.2 and iqr > 0
            else "High zero fraction or low IQR → near detection floor, rank correlation compressed"),
    }


def diag_power(name, n):
    """D2: Power analysis — α=0.05 two-sided minimum detectable |ρ|"""
    df = n - 2
    # t critical value: inverting via betai is inconvenient; use approximation + exact lookup table combination
    # Small samples use known critical value table (df<=30) + large-sample normal approximation
    tcrit_table = {
        22: 2.074, 23: 2.069, 30: 2.042,  # n=24→df22, n=32→df30
    }
    if df in tcrit_table:
        tcrit = tcrit_table[df]
    else:
        # Large-sample approximation (df>30): t ≈ z_{0.975} = 1.96
        tcrit = 1.96
    # |r| = t / sqrt(t^2 + df)
    min_r = tcrit / math.sqrt(tcrit * tcrit + df)
    return {
        "Cohort": name, "Diag": "D2_power", "Metric": f"n={n}",
        "Value": f"alpha=0.05,two-sided,min_detectable_|rho|={min_r:.3f}",
        "Interpretation": (
            f"At this sample size, can detect at most |rho|≥{min_r:.2f}; true effects 0.2–0.3 likely fall into null interval"),
    }


def diag_scale(name, zp3, score):
    """D3: expression scale sensitivity"""
    out = []
    # original
    r_raw, _ = spearman(zp3, score)
    # log1p
    zp3_log = [math.log1p(max(0.0, v)) for v in zp3]
    r_log, _ = spearman(zp3_log, score)
    # binary (detected / not detected)
    zp3_bin = [1.0 if v and v > 0 else 0.0 for v in zp3]
    r_bin, _ = spearman(zp3_bin, score)
    # high-expression stratification (≥Q3 vs <Q3)
    q3 = quantile(zp3, 0.75)
    zp3_hi = [1.0 if v >= q3 else 0.0 for v in zp3]
    r_hi, _ = spearman(zp3_hi, score)
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "raw",
        "Value": f"rho={r_raw:+.3f}", "Interpretation": "raw expression",
    })
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "log1p",
        "Value": f"rho={r_log:+.3f}", "Interpretation": "log1p transformation",
    })
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "detect_binary",
        "Value": f"rho={r_bin:+.3f}", "Interpretation": "detected/not detected binary",
    })
    out.append({
        "Cohort": name, "Diag": "D3_scale_M2", "Metric": "high_vs_low",
        "Value": f"rho={r_hi:+.3f}", "Interpretation": "high expression (≥Q3) vs low expression",
    })
    return out


def diag_jackknife(name, zp3, score):
    """D4: jackknife stability — M2 ρ fluctuation after removing each sample one by one"""
    n = len(zp3)
    rhos = []
    for i in range(n):
        x = zp3[:i] + zp3[i + 1:]
        y = score[:i] + score[i + 1:]
        r, _ = spearman(x, y)
        rhos.append(r)
    rmin, rmax = min(rhos), max(rhos)
    rmean = sum(rhos) / n
    # compare with full sample
    r_full, _ = spearman(zp3, score)
    spread = rmax - rmin
    all_neg = rmin < 0 and rmax < 0
    return {
        "Cohort": name, "Diag": "D4_jackknife_M2", "Metric": f"n={n}",
        "Value": f"full_rho={r_full:+.3f},leave1out_range=[{rmin:+.3f},{rmax:+.3f}],"
                 f"spread={spread:.3f},mean={rmean:+.3f}",
        "Interpretation": (
            "Direction consistently negative and moderate magnitude sensitivity → not flipped by individual samples"
            if spread < 0.25 and all_neg
            else "Wide range → result sensitive to individual samples"),
    }


def main():
    assert os.path.isfile(GSE77530), f"Missing: {GSE77530}"
    assert os.path.isfile(GSE113474), f"Missing: {GSE113474}"

    cohorts = [
        ("GSE77530_MDAnderson", GSE77530, "\t", "RPKM"),
        ("GSE113474_NYU", GSE113474, ",", "norm_counts"),
    ]
    rows_out = []
    print("=== A3 Freeze⑦ External Null Diagnostics ===")
    for name, path, delim, unit in cohorts:
        samples, mat = load_matrix(path, delim)
        if "ZP3" not in mat:
            print(f"  {name}: ZP3 missing, skipping"); continue
        zp3 = mat["ZP3"]
        n = len(samples)
        print(f"\n{name}: n={n}, unit={unit}")

        # D1
        d1 = diag_detection(name, zp3)
        rows_out.append(d1)
        print(f"  D1 detection quality: {d1['Value']} → {d1['Interpretation']}")

        # D2
        d2 = diag_power(name, n)
        rows_out.append(d2)
        print(f"  D2 power: {d2['Value']}")

        # M2 immune score (baseline for diagnostic scale + jackknife)
        m2_genes = IMMUNE['M2_Macrophage']
        score, used = zscore_consensus(mat, m2_genes)
        if score is None:
            print(f"  M2 scoring genes insufficient, skipping D3/D4"); continue

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
    print(f"\nFrozen table: {OUT_CSV}")

    # ---- Self-check: diagnostics do not affect main conclusions, only verify format and reasonableness ----
    print("\n=== Self-check ===")
    ok = True
    # Self-check 1: all rows have interpretation
    for r in rows_out:
        if not r.get("Interpretation"):
            print(f"FAIL: missing interpretation {r}"); ok = False
    # Self-check 2: D2 minimum detectable |rho| in reasonable range (0.3–0.6)
    for r in rows_out:
        if r["Diag"] == "D2_power":
            s = r["Value"]
            lo = s.rfind("=") + 1
            mr = float(s[lo:].strip())
            if not (0.3 <= mr <= 0.6):
                print(f"FAIL: D2 minimum |rho| abnormal {mr}"); ok = False
    # Self-check 3: D1 zero proportion non-negative
    for r in rows_out:
        if r["Diag"] == "D1_detection":
            if "zero=" not in r["Value"]:
                print(f"FAIL: D1 format abnormal"); ok = False
    print("  Note: this script only does diagnostics, does not recompute the main null; main null see a3_external_gbm.csv")
    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
