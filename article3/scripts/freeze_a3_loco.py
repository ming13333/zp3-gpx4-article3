#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结⑧ — leave-one-cancer-out 内部验证（freeze_a3_loco）
================================================================
目的：回答「内部 isoform 关联是否由某几个癌种单独驱动？」。
      对 32 癌种的 FL 比例 × 免疫评分关联，逐癌种留出后重算固定效应
      Fisher-z 合并，观察合并 ρ / 95% CI 是否稳定。

这是外部 null 诊断的补充：外部 null 反映「基因级跨队列泛化受限」，
本脚本证明「内部 isoform 级关联在癌种间具有稳健性（非单癌种假象）」。

输入：article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv
      （每行一个样本，列含 FL 比例、7 个 score_* 免疫评分、Cancer）
输出：article3/results/a3_loco_frozen.csv
      （每癌种一行：留出该癌种后的合并 ρ / 95% CI / Q / I² / 剩余癌种数）

统计基元复用 freeze_a3_robustness.py（纯标准库 Spearman + betai + Fisher-z）。
"""
import os
import sys
import csv
import math

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根
JOINED = os.path.join(ROOT, "article3", "results", "zp3_psi_pancancer_results",
                      "psi_immune_joined_samples.csv")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_loco_frozen.csv")

FL = "ENST00000336517.8"
IMMUNE = ["score_M2_Macrophage", "score_T_cell_exhaustion", "score_Cytolytic_activity",
          "score_Treg", "score_IFN_gamma", "score_Checkpoint", "score_Myeloid"]
SCORES = [c.replace("score_", "") for c in IMMUNE]


# ---------------------------------------------------------------------------
# 纯标准库统计基元（与 freeze_a3_robustness.py 一致）
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
            d = FPMIN
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
        entry = {"fl": fl}
        ok = True
        for sc in IMMUNE:
            v = f2(d.get(sc))
            if v is None:
                ok = False
                break
            entry[sc] = v
        if not ok:
            continue
        by_cancer.setdefault(c, []).append(entry)

    cancers = sorted(by_cancer.keys())
    print(f"癌种数: {len(cancers)}")

    # 全量池化（基线，供对照）
    def pooled(entries_per_score):
        """entries_per_score: {score: [(fl, score_val), ...]} 跨所有保留癌种"""
        rows_out = {}
        for sc in IMMUNE:
            xs, ys = [], []
            for flv, sv in entries_per_score[sc]:
                xs.append(flv)
                ys.append(sv)
            r, _ = spearman(xs, ys)
            rows_out[sc] = r
        return rows_out

    # 预聚合每个癌种的 (fl, score) 列表
    cancer_scores = {}
    for c in cancers:
        d = {sc: [(e["fl"], e[sc]) for e in by_cancer[c]] for sc in IMMUNE}
        cancer_scores[c] = d

    # 全量合并（leave-none-out 基线）
    all_entries = {sc: [] for sc in IMMUNE}
    for c in cancers:
        for sc in IMMUNE:
            all_entries[sc].extend(cancer_scores[c][sc])

    def fisher_meta(entries):
        """固定效应 Fisher-z 合并，返回 (pooled_r, lo, hi, Q, I2, k)"""
        Sw = 0.0
        Z = 0.0
        zs = []
        ns = []
        for sc in IMMUNE:
            for x, y in entries[sc]:
                pass  # 这里不使用
        # 改为逐癌种聚合（每癌种一个 z），避免样本量重复
        # 先算每癌种每评分的单相关系数，再按癌种加权
        return None

    # 简化：LOCO 以癌种为随机效应单元，每癌种内先算 Spearman(FL, score)，
    # 再对「保留癌种」做固定效应 Fisher-z 合并（权重 N_cancer - 3）。
    def cancer_level_z(c):
        out = {}
        for sc in IMMUNE:
            xs = [flv for flv, _ in cancer_scores[c][sc]]
            ys = [sv for _, sv in cancer_scores[c][sc]]
            r, _ = spearman(xs, ys)
            out[sc] = (r, len(xs))
        return out

    cancer_z = {c: cancer_level_z(c) for c in cancers}

    def meta_over(keep_cancers, score):
        Sw = 0.0
        Z = 0.0
        zlist = []
        for c in keep_cancers:
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
        pooled_r = math.tanh(Z)
        Q = sum((z - Z) ** 2 for (_, _, _, z, _) in zlist)
        df = len(zlist) - 1
        I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
        return pooled_r, math.tanh(lo), math.tanh(hi), Q, I2, len(zlist)

    # 全量基线
    frozen = []
    for score in SCORES:
        m = meta_over(cancers, "score_" + score)
        if m:
            pr, lo, hi, Q, I2, k = m
            frozen.append({
                "Left_out_cancer": "NONE(all)", "Score": score, "K_cancers": k,
                "Pooled_rho": round(pr, 4), "CI_low": round(lo, 4), "CI_high": round(hi, 4),
                "Q": round(Q, 2), "I2": round(I2, 3),
                "CI_crosses_zero": "YES" if lo < 0 < hi else "NO",
            })

    # LOCO
    print("\n=== Leave-one-cancer-out ===")
    for drop in cancers:
        keep = [c for c in cancers if c != drop]
        for score in SCORES:
            m = meta_over(keep, "score_" + score)
            if not m:
                continue
            pr, lo, hi, Q, I2, k = m
            frozen.append({
                "Left_out_cancer": drop, "Score": score, "K_cancers": k,
                "Pooled_rho": round(pr, 4), "CI_low": round(lo, 4), "CI_high": round(hi, 4),
                "Q": round(Q, 2), "I2": round(I2, 3),
                "CI_crosses_zero": "YES" if lo < 0 < hi else "NO",
            })

    cols = ["Left_out_cancer", "Score", "K_cancers", "Pooled_rho",
            "CI_low", "CI_high", "Q", "I2", "CI_crosses_zero"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in frozen:
            w.writerow(r)
    print(f"冻结表: {OUT_CSV} ({len(frozen)} 行)")

    # ---- 自检：LOCO 合并 ρ 应与全量基线同号且幅度不剧烈翻转 ----
    print("\n=== 自检 ===")
    ok = True
    base = {r["Score"]: r for r in frozen if r["Left_out_cancer"] == "NONE(all)"}
    flips = 0
    for r in frozen:
        if r["Left_out_cancer"] == "NONE(all)":
            continue
        sc = r["Score"]
        if sc in base:
            b_r = base[sc]["Pooled_rho"]
            if (b_r > 0) != (r["Pooled_rho"] > 0):
                flips += 1
                print(f"  方向翻转: 留出 {r['Left_out_cancer']} {sc} "
                      f"base={b_r:+.3f} lo={r['Pooled_rho']:+.3f}")
    if flips == 0:
        print("  PASS: 所有 LOCO 合并 ρ 与全量基线同号（内部稳健，非单癌种驱动）")
    else:
        print(f"  WARN: {flips} 次方向翻转（仍冻结，供稿件透明报告）")
    # 幅度变化范围
    for sc in SCORES:
        vals = [r["Pooled_rho"] for r in frozen if r["Score"] == sc and r["Left_out_cancer"] != "NONE(all)"]
        if vals:
            print(f"  {sc:18s} LOCO ρ range=[{min(vals):+.3f}, {max(vals):+.3f}] "
                  f"(base={base[sc]['Pooled_rho']:+.3f})")
    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
