#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 Freeze ⑨ — Internal cross-cancer transportability validation (freeze_a3_transportability)
================================================================
Objective: Without relying on external permission-controlled data, solidify the "internal generalizability" validation:
  A. LOCO (leave-one-cancer-out) — already frozen in freeze_a3_loco.py, recomputed here and provided
     comparison of direction/magnitude with the full baseline (M2/Myeloid two headline scores).
  B. L2CO (leave-two-cancer-out) — leave out any two cancer types simultaneously, test whether the combined ρ remains robust
     (a more stringent sensitivity analysis than LOCO, C(32,2)=496 pairs).
  C. Preset held-out split — fix seed=42, randomly split train/val 70/30 by cancer type (split by cancer type
     to avoid leakage of same-cancer-type samples across sets); train set estimates, val set performs independent direction test.
     Note: the split is preset (seed fixed, not selected based on results); it is internal transportability
     rather than external independent validation (the external gene-level null has been frozen at a3_external_gbm.csv).

Only run for the two headline scores M2_Macrophage and Myeloid (consistent with the robustness meta definition).

Input: article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv
# Output: article3/results/a3_transportability_frozen.csv
# Statistical primitives: pure standard library (same as freeze_a3_loco.py / freeze_a3_robustness.py).
"""
import os
import sys
import csv
import math
import random
import itertools

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 levels = project root
JOINED = os.path.join(ROOT, "article3", "results", "zp3_psi_pancancer_results",
                      "psi_immune_joined_samples.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_transportability_frozen.csv")

FL = "ENST00000336517.8"
SCORES = ["M2_Macrophage", "Myeloid"]
IMMUNE = ["score_" + s for s in SCORES]


# ---------------------------------------------------------------------------
# Pure standard library statistical primitives (consistent implementation)
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


def load():
    with open(JOINED, newline="") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
def main():
    assert os.path.isfile(JOINED), f"Input missing: {JOINED}"
    rows = load()
    print(f"Joined samples: {len(rows)}")

    # Group by cancer type
    by_cancer = {}
    for d in rows:
        c = d.get("Cancer")
        if not c:
            continue
        fl = f2(d.get(FL))
        if fl is None:
            continue
        ok = True
        vals = {"fl": fl}
        for sc in IMMUNE:
            v = f2(d.get(sc))
            if v is None:
                ok = False
                break
            vals[sc] = v
        if not ok:
            continue
        by_cancer.setdefault(c, []).append(vals)

    cancers = sorted(by_cancer.keys())
    print(f"Number of cancer types: {len(cancers)}")

    # Spearman(FL, score) for each cancer type and score
    cancer_z = {}
    for c in cancers:
        cancer_z[c] = {}
        for sc in IMMUNE:
            xs = [e["fl"] for e in by_cancer[c]]
            ys = [e[sc] for e in by_cancer[c]]
            r, _ = spearman(xs, ys)
            cancer_z[c][sc] = (r, len(xs))

    def meta_over(keep, score):
        Sw, Z = 0.0, 0.0
        zlist = []
        for c in keep:
            r, n = cancer_z[c][score]
            if n < 4:
                continue
            w = n - 3
            z = math.atanh(max(-0.9999, min(0.9999, r)))
            Sw += w
            Z += w * z
            zlist.append((c, r, n, z, w))
        if Sw == 0:
            return None
        Z /= Sw
        se = 1.0 / math.sqrt(Sw)
        lo, hi = Z - 1.96 * se, Z + 1.96 * se
        pr = math.tanh(Z)
        Q = sum(w * (z - Z) ** 2 for (_, _, _, z, w) in zlist)
        df = len(zlist) - 1
        I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
        return pr, math.tanh(lo), math.tanh(hi), Q, I2, len(zlist)

    frozen = []

    def add_row(analysis, score, desc, k, pr, lo, hi, Q, I2,
                train_r=None, val_r=None):
        frozen.append({
            "Analysis": analysis, "Score": score, "Split": desc, "K": k,
            "Pooled_rho": round(pr, 4), "CI_low": round(lo, 4), "CI_high": round(hi, 4),
            "Q": round(Q, 2), "I2": round(I2, 3),
            "Train_rho": "" if train_r is None else round(train_r, 4),
            "Val_rho": "" if val_r is None else round(val_r, 4),
            "Same_direction": "" if train_r is None else ("YES" if (train_r >= 0) == (val_r >= 0) else "NO"),
            "Val_CI_crosses_zero": "" if val_r is None else ("YES" if lo < 0 < hi else "NO"),
        })

    # ---- Full baseline ----
    print("\n=== Full baseline ===")
    for sc in IMMUNE:
        m = meta_over(cancers, sc)
        if m:
            pr, lo, hi, Q, I2, k = m
            add_row("BASELINE", sc, "all_32", k, pr, lo, hi, Q, I2)
            print(f"  {sc:16s} rho={pr:+.3f} (CI {lo:+.3f}~{hi:+.3f}) I²={I2*100:.0f}%")

    # ---- A. LOCO ----
    print("\n=== A. LOCO (leave-one-cancer-out) ===")
    for drop in cancers:
        keep = [c for c in cancers if c != drop]
        for sc in IMMUNE:
            m = meta_over(keep, sc)
            if m:
                pr, lo, hi, Q, I2, k = m
                add_row("LOCO", sc, f"leave:{drop}", k, pr, lo, hi, Q, I2)
    print(f"  LOCO row count: {sum(1 for r in frozen if r['Analysis']=='LOCO')}")

    # ---- B. L2CO ----
    print("\n=== B. L2CO (leave-two-cancer-out) ===")
    pairs = list(itertools.combinations(cancers, 2))
    for c1, c2 in pairs:
        keep = [c for c in cancers if c != c1 and c != c2]
        for sc in IMMUNE:
            m = meta_over(keep, sc)
            if m:
                pr, lo, hi, Q, I2, k = m
                add_row("L2CO", sc, f"leave:{c1}+{c2}", k, pr, lo, hi, Q, I2)
    print(f"  L2CO rows: {sum(1 for r in frozen if r['Analysis']=='L2CO')} "
          f"({len(pairs)} pairs × {len(IMMUNE)} scores)")

    # ---- C. Preset held-out split (seed=42, by cancer type 70/30) ----
    print("\n=== C. HELDOUT (prespecified, seed=42, by cancer type 70/30) ===")
    rng = random.Random(42)
    shuffled = list(cancers)
    rng.shuffle(shuffled)
    n_train = int(round(len(shuffled) * 0.7))
    train_set, val_set = set(shuffled[:n_train]), set(shuffled[n_train:])
    print(f"  train {len(train_set)} cancer types / val {len(val_set)} cancer types")
    print(f"  train: {sorted(train_set)}")
    print(f"  val:   {sorted(val_set)}")
    for sc in IMMUNE:
        m_tr = meta_over(sorted(train_set), sc)
        m_va = meta_over(sorted(val_set), sc)
        if m_tr and m_va:
            pr_t, lo_t, hi_t, Q_t, I2_t, k_t = m_tr
            pr_v, lo_v, hi_v, Q_v, I2_v, k_v = m_va
            add_row("HELDOUT", sc, "train32-70_seed42", k_t, pr_t, lo_t, hi_t, Q_t, I2_t,
                    train_r=pr_t, val_r=pr_v)
            # val record its CI on a separate line
            frozen.append({
                "Analysis": "HELDOUT_VAL", "Score": sc, "Split": "val30_seed42", "K": k_v,
                "Pooled_rho": round(pr_v, 4), "CI_low": round(lo_v, 4), "CI_high": round(hi_v, 4),
                "Q": round(Q_v, 2), "I2": round(I2_v, 3),
                "Train_rho": round(pr_t, 4), "Val_rho": round(pr_v, 4),
                "Same_direction": "YES" if (pr_t >= 0) == (pr_v >= 0) else "NO",
                "Val_CI_crosses_zero": "YES" if lo_v < 0 < hi_v else "NO",
            })
            print(f"  {sc:16s} train rho={pr_t:+.3f} val rho={pr_v:+.3f} "
                  f"same_dir={'YES' if (pr_t>=0)==(pr_v>=0) else 'NO'} "
                  f"val_CI_cross_zero={'YES' if lo_v<0<hi_v else 'NO'}")

    cols = ["Analysis", "Score", "Split", "K", "Pooled_rho", "CI_low", "CI_high",
            "Q", "I2", "Train_rho", "Val_rho", "Same_direction", "Val_CI_crosses_zero"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in frozen:
            w.writerow(r)
    print(f"\nFrozen table: {OUT_CSV} ({len(frozen)} rows)")

    # ---- self-check ----
    print("\n=== Self-check ===")
    ok = True
    base = {r["Score"]: r for r in frozen if r["Analysis"] == "BASELINE"}
    for sc in IMMUNE:
        # L2CO should not have direction reversal (relative to baseline)
        flips_l2co = 0
        for r in frozen:
            if r["Analysis"] == "L2CO" and r["Score"] == sc:
                b = base[sc]["Pooled_rho"]
                if (b >= 0) != (r["Pooled_rho"] >= 0):
                    flips_l2co += 1
        print(f"  L2CO {sc:16s} direction flips {flips_l2co}/{len(pairs)} pairs")
        # HELDOUT direction consistent + val CI crossing 0 is expected (small sample, low power, honest report)
        ho = [r for r in frozen if r["Analysis"] == "HELDOUT" and r["Score"] == sc]
        if ho and ho[0]["Same_direction"] != "YES":
            print(f"  FAIL: HELDOUT {sc} train/val direction mismatch"); ok = False
    print("  Note: val set has only ~9 cancer types, wide CI is expected; HELDOUT is internal transportability, "
          "not external independent validation (the latter see a3_external_gbm.csv)")
    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
