#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结⑤ — 组成数据与稳健性分析（freeze_a3_robustness）
================================================================
目的：正面回应审稿人对「组成数据 / compositional data」与「样本级验证薄弱」
两项核心质疑，从已冻结的输入独立复算，产出两张冻结表：
        - a3_robustness_frozen.csv  （比例 metric 敏感性 + 组成耦合）
        - a3_robustness_meta.csv     （癌种分层 Fisher-z 荟萃 + 异质性）

四大分析：
  1) Baseline：FL 比例 vs 7 个免疫评分（Spearman/Pearson）—— 与 mixed-model 方向对照
  2) Compositional contrast：log(FL/RI) 配对 log-ratio 替代比例，验证方向稳健性
     （log-ratio 是组成数据最轻量的 ILR 控制，规避单纯比例的组成约束误导）
  3) FL–RI coupling：Spearman(FL, RI) 量化两主导 isoform 比例的负相关耦合
     （说明比例是组成约束下的相对量，不能直接解读为「独立丰度」）
  4) Low-signal filter：仅保留主导 isoform 比例 ≥ 0.5 的样本（可信定量代理），
     复算 FL 比例关联，检验低信号样本是否驱动结论

荟萃（meta）：
  对 headline 免疫评分（M2_Macrophage、Myeloid），按癌种分层计算 Spearman(FL, score)
  （每癌种 N≥30），固定效应 Fisher-z 合并：Z=Σw_i z_i / Σw_i，w_i=N_i−3；
  报告合并 ρ、95% CI、Cochran Q、I²。

实现要点（纯标准库，无第三方依赖，最大化可复现性）：
  - Spearman = 对平均秩（处理并列）做 Pearson
  - 相关 p 值：大样本 t 近似经正则不完全 Beta（betai）精确校正
  - 比例极小值加 eps 防 log 下溢

输入：
  - article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv
输出：
  - article3/results/a3_robustness_frozen.csv
  - article3/results/a3_robustness_meta.csv
"""
import os
import sys
import csv
import math

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根（实测验证）
JOINED = os.path.join(ROOT, "article3", "results", "zp3_psi_pancancer_results",
                      "psi_immune_joined_samples.csv")
OUT1 = os.path.join(ROOT, "article3", "results", "a3_robustness_frozen.csv")
OUT2 = os.path.join(ROOT, "article3", "results", "a3_robustness_meta.csv")

assert os.path.isfile(JOINED), f"输入缺失: {JOINED}"

FL = "ENST00000336517.8"
RI = "ENST00000466960.5"
ISO_COLS = [FL, RI, "ENST00000394860.3", "ENST00000467555.1",
            "ENST00000394857.7", "ENST00000416245.5", "ENST00000479793.5"]
IMMUNE = ["score_M2_Macrophage", "score_T_cell_exhaustion", "score_Cytolytic_activity",
          "score_Treg", "score_IFN_gamma", "score_Checkpoint", "score_Myeloid"]
EPS = 1e-9


# ---------------------------------------------------------------------------
# 统计基元（纯标准库）
# ---------------------------------------------------------------------------
def f2(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rankdata(vals):
    """平均秩（处理并列）。vals: list[float]"""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # 1-based 平均秩
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
    r = sxy / math.sqrt(sxx * syy)
    if r > 1.0:
        r = 1.0
    if r < -1.0:
        r = -1.0
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
    """返回 (r, p_two_sided)，基于 t 分布精确 p 值。"""
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
    r, p = corr_test(rx, ry)
    return r, p


def label(s):
    return s.replace("score_", "")


# ---------------------------------------------------------------------------
def load():
    rows = []
    with open(JOINED, newline="") as f:
        rd = csv.DictReader(f)
        for d in rd:
            rows.append(d)
    return rows


def paired(rows, key_x, key_y):
    xs, ys = [], []
    for d in rows:
        a = f2(d.get(key_x))
        b = f2(d.get(key_y))
        if a is None or b is None:
            continue
        xs.append(a)
        ys.append(b)
    return xs, ys


def main():
    rows = load()
    print(f"Joined 样本: {len(rows)}")

    frozen = []  # tidy long

    # ---- (1) Baseline: FL 比例 vs 7 免疫评分 ----
    print("\n=== (1) Baseline: FL 比例 vs 免疫评分 ===")
    baseline = {}
    for s in IMMUNE:
        x, y = paired(rows, FL, s)
        r, p = spearman(x, y)
        rp, pp = corr_test(x, y)
        baseline[s] = (r, p)
        frozen.append({
            "Analysis": "FL_proportion_baseline", "Metric": "proportion",
            "Immune_score": label(s), "N": len(x),
            "Spearman_rho": round(r, 4), "Spearman_p": f"{p:.3e}",
            "Pearson_r": round(rp, 4), "Pearson_p": f"{pp:.3e}", "Note": "vs FL canonical proportion",
        })
        print(f"  {label(s):22s} n={len(x):5d}  rho={r:+.3f} (p={p:.1e})  r={rp:+.3f}")

    # ---- (2) Compositional contrast: log(FL/RI) ----
    print("\n=== (2) Compositional contrast: log(FL/RI) ===")
    logratio = {}
    for s in IMMUNE:
        xs, ys = [], []
        for d in rows:
            a = f2(d.get(FL)); b = f2(d.get(RI)); c = f2(d.get(s))
            if a is None or b is None or c is None:
                continue
            lr = math.log((a + EPS) / (b + EPS))
            xs.append(lr)
            ys.append(c)
        r, p = spearman(xs, ys)
        rp, pp = corr_test(xs, ys)
        logratio[s] = (r, p)
        frozen.append({
            "Analysis": "log_FL_RI_ratio", "Metric": "compositional_contrast",
            "Immune_score": label(s), "N": len(xs),
            "Spearman_rho": round(r, 4), "Spearman_p": f"{p:.3e}",
            "Pearson_r": round(rp, 4), "Pearson_p": f"{pp:.3e}", "Note": "pairwise log-ratio FL/RI",
        })
        print(f"  {label(s):22s} n={len(xs):5d}  rho={r:+.3f} (p={p:.1e})  r={rp:+.3f}")

    # ---- (3) FL–RI coupling ----
    print("\n=== (3) FL–RI compositional coupling ===")
    fx, rx = paired(rows, FL, RI)
    r_coup, p_coup = spearman(fx, rx)
    frozen.append({
        "Analysis": "FL_RI_coupling", "Metric": "compositional",
        "Immune_score": "FL_vs_RI", "N": len(fx),
        "Spearman_rho": round(r_coup, 4), "Spearman_p": f"{p_coup:.3e}",
        "Pearson_r": "", "Pearson_p": "", "Note": "component coupling within compositional sum",
    })
    print(f"  Spearman(FL, RI) = {r_coup:+.3f} (p={p_coup:.1e})")

    # ---- (4) Low-signal filter: 主导 isoform 比例 ≥ 0.5 ----
    print("\n=== (4) Low-signal filter (max isoform proportion >= 0.5) ===")
    kept = []
    for d in rows:
        vals = [f2(d.get(c)) for c in ISO_COLS]
        if any(v is None for v in vals):
            continue
        if max(vals) >= 0.5:
            kept.append(d)
    print(f"  保留样本: {len(kept)} / {len(rows)}")
    for s in IMMUNE:
        x, y = paired(kept, FL, s)
        r, p = spearman(x, y)
        rp, pp = corr_test(x, y)
        frozen.append({
            "Analysis": "LowSignal_filter_maxIso_ge_0.5", "Metric": "proportion_filtered",
            "Immune_score": label(s), "N": len(x),
            "Spearman_rho": round(r, 4), "Spearman_p": f"{p:.3e}",
            "Pearson_r": round(rp, 4), "Pearson_p": f"{pp:.3e}",
            "Note": f"filtered from {len(rows)} to {len(kept)}",
        })
        print(f"  {label(s):22s} n={len(x):5d}  rho={r:+.3f} (p={p:.1e})")

    # ---- 写入 frozen CSV ----
    cols = ["Analysis", "Metric", "Immune_score", "N", "Spearman_rho",
            "Spearman_p", "Pearson_r", "Pearson_p", "Note"]
    with open(OUT1, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in frozen:
            w.writerow(r)
    print(f"\n冻结表1: {OUT1}")

    # ---- Meta-analysis: 按癌种分层 Fisher-z 合并 ----
    print("\n=== Meta-analysis (cancer-stratified Fisher-z) ===")
    meta_cols = ["Score", "Cancer", "N", "Spearman_rho", "Spearman_p",
                 "Fisher_z", "Weight", "Q", "I2"]
    meta_rows = []
    for s in ["score_M2_Macrophage", "score_Myeloid"]:
        by_cancer = {}
        for d in rows:
            c = d.get("Cancer")
            a = f2(d.get(FL)); b = f2(d.get(s))
            if a is None or b is None or not c:
                continue
            by_cancer.setdefault(c, {"x": [], "y": []})
            by_cancer[c]["x"].append(a)
            by_cancer[c]["y"].append(b)
        per = []
        for c, dd in by_cancer.items():
            if len(dd["x"]) < 30:
                continue
            r, p = spearman(dd["x"], dd["y"])
            z = math.atanh(max(-0.9999, min(0.9999, r)))
            per.append((c, len(dd["x"]), r, p, z, len(dd["x"]) - 3))
        if len(per) < 2:
            print(f"  {label(s)}: 可合并癌种 <2，跳过")
            continue
        Sw = sum(w for *_, w in per)
        Z = sum(z * w for *_, z, w in per) / Sw
        se = 1.0 / math.sqrt(Sw)
        lo, hi = Z - 1.96 * se, Z + 1.96 * se
        pooled_r = math.tanh(Z)
        ci_lo, ci_hi = math.tanh(lo), math.tanh(hi)
        Q = sum(w * (z - Z) ** 2 for *_, z, w in per)
        df = len(per) - 1
        I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
        # pooled p (normal approx on Z)
        pooled_p = 2.0 * (0.5 * math.erfc(abs(Z) / math.sqrt(2)))
        for c, n, r, p, z, w in per:
            meta_rows.append({
                "Score": label(s), "Cancer": c, "N": n, "Spearman_rho": round(r, 4),
                "Spearman_p": f"{p:.3e}", "Fisher_z": round(z, 4), "Weight": w,
                "Q": "", "I2": "",
            })
        meta_rows.append({
            "Score": label(s), "Cancer": "POOLED_fixed_effect", "N": int(Sw + 3 * len(per)),
            "Spearman_rho": round(pooled_r, 4), "Spearman_p": f"{pooled_p:.3e}",
            "Fisher_z": round(Z, 4), "Weight": round(Sw, 1),
            "Q": round(Q, 2), "I2": round(I2, 3),
        })
        print(f"  {label(s):12s} k={len(per)} pooled ρ={pooled_r:+.3f} "
              f"(95% CI {ci_lo:+.3f}–{ci_hi:+.3f}) Q={Q:.1f} I²={I2*100:.0f}%")
    with open(OUT2, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=meta_cols)
        w.writeheader()
        for r in meta_rows:
            w.writerow(r)
    print(f"冻结表2: {OUT2}")

    # ---- 自检 ----
    print("\n=== 自评 ===")
    ok = True
    b_m2 = baseline["score_M2_Macrophage"][0]
    l_m2 = logratio["score_M2_Macrophage"][0]
    if b_m2 <= 0:
        ok = False; print("  FAIL: baseline FL vs M2 应为正")
    else:
        print(f"  PASS: baseline FL vs M2 ρ={b_m2:+.3f} (与 mixed-model β=+0.28 方向一致)")
    if l_m2 <= 0:
        ok = False; print("  FAIL: log(FL/RI) vs M2 应为正")
    else:
        print(f"  PASS: log(FL/RI) vs M2 ρ={l_m2:+.3f} (组成控制后方向不变)")
    if r_coup >= 0:
        ok = False; print("  FAIL: FL-RI 应负相关（组成耦合）")
    else:
        print(f"  PASS: FL-RI Spearman={r_coup:+.3f} (组成耦合，比例非独立丰度)")
    # meta M2 pooled
    m2_pooled = [r for r in meta_rows if r["Score"] == "M2_Macrophage" and r["Cancer"] == "POOLED_fixed_effect"]
    if m2_pooled:
        pr = float(m2_pooled[0]["Spearman_rho"])
        if pr <= 0:
            ok = False; print("  FAIL: M2 荟萃合并 ρ 应为正")
        else:
            print(f"  PASS: M2 跨癌种荟萃 ρ={pr:+.3f}")

    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
