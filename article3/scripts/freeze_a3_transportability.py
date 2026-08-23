#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结⑨ — 内部跨癌种迁移性验证（freeze_a3_transportability）
================================================================
目的：在不依赖外部权限数据的前提下，把「内部泛化性」验证做扎实：
  A. LOCO（leave-one-cancer-out）—— 已在 freeze_a3_loco.py 冻结，此处重算并给出
     与全量基线的方向/幅度对照（M2/Myeloid 两个 headline 评分）。
  B. L2CO（leave-two-cancer-out）—— 任意两个癌种同时留出，检验合并 ρ 是否仍稳健
     （比 LOCO 更严苛的敏感性分析，C(32,2)=496 对）。
  C. 预设 held-out split —— 固定 seed=42、按癌种 70/30 随机分 train/val（按癌种分割
     避免同癌种样本跨集泄漏）；train 集估计，val 集做独立方向检验。
     强调：split 为预设（seed 固定、不按结果挑选），是内部 transportability
     而非外部独立验证（外部 gene-level null 已冻结在 a3_external_gbm.csv）。

仅对 M2_Macrophage 与 Myeloid 两个 headline 评分执行（与 robustness meta 口径一致）。

输入：article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv
输出：article3/results/a3_transportability_frozen.csv
统计基元：纯标准库（同 freeze_a3_loco.py / freeze_a3_robustness.py）。
"""
import os
import sys
import csv
import math
import random
import itertools

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根
JOINED = os.path.join(ROOT, "article3", "results", "zp3_psi_pancancer_results",
                      "psi_immune_joined_samples.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_transportability_frozen.csv")

FL = "ENST00000336517.8"
SCORES = ["M2_Macrophage", "Myeloid"]
IMMUNE = ["score_" + s for s in SCORES]


# ---------------------------------------------------------------------------
# 纯标准库统计基元（一致实现）
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
    assert os.path.isfile(JOINED), f"输入缺失: {JOINED}"
    rows = load()
    print(f"Joined 样本: {len(rows)}")

    # 按癌种分组
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
    print(f"癌种数: {len(cancers)}")

    # 每癌种每评分 Spearman(FL, score)
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

    # ---- 全量基线 ----
    print("\n=== 全量基线 ===")
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
    print(f"  LOCO 行数: {sum(1 for r in frozen if r['Analysis']=='LOCO')}")

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
    print(f"  L2CO 行数: {sum(1 for r in frozen if r['Analysis']=='L2CO')} "
          f"({len(pairs)} 对 × {len(IMMUNE)} 评分)")

    # ---- C. 预设 held-out split（seed=42，按癌种 70/30）----
    print("\n=== C. HELDOUT (prespecified, seed=42, 按癌种 70/30) ===")
    rng = random.Random(42)
    shuffled = list(cancers)
    rng.shuffle(shuffled)
    n_train = int(round(len(shuffled) * 0.7))
    train_set, val_set = set(shuffled[:n_train]), set(shuffled[n_train:])
    print(f"  train {len(train_set)} 癌种 / val {len(val_set)} 癌种")
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
            # val 单独一行记录其 CI
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
    print(f"\n冻结表: {OUT_CSV} ({len(frozen)} 行)")

    # ---- 自检 ----
    print("\n=== 自检 ===")
    ok = True
    base = {r["Score"]: r for r in frozen if r["Analysis"] == "BASELINE"}
    for sc in IMMUNE:
        # L2CO 不应出现方向翻转（相对 baseline）
        flips_l2co = 0
        for r in frozen:
            if r["Analysis"] == "L2CO" and r["Score"] == sc:
                b = base[sc]["Pooled_rho"]
                if (b >= 0) != (r["Pooled_rho"] >= 0):
                    flips_l2co += 1
        print(f"  L2CO {sc:16s} 方向翻转 {flips_l2co}/{len(pairs)} 对")
        # HELDOUT 方向一致 + val CI 跨 0 为预期（样本少功效低，诚实报告）
        ho = [r for r in frozen if r["Analysis"] == "HELDOUT" and r["Score"] == sc]
        if ho and ho[0]["Same_direction"] != "YES":
            print(f"  FAIL: HELDOUT {sc} train/val 方向不一致"); ok = False
    print("  注：val 集仅 ~9 癌种，CI 较宽属预期；HELDOUT 为内部 transportability，"
          "非外部独立验证（后者见 a3_external_gbm.csv）")
    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()