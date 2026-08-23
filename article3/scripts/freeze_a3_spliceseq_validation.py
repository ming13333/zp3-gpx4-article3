#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结② — SpliceSeq 代理验证（freeze_a3_spliceseq_validation）
================================================================
目的：从独立输入（SpliceSeq PSI 下载 + 比例矩阵）复算 A3 图 2 全部统计量，
      产出两张冻结表：
        - a3_spliceseq_ecological.csv  （癌种级生态学 ρ/P + LOO + bootstrap CI）
        - a3_spliceseq_samplelevel.csv （GBM/LGG 样本级各 AP 事件 Spearman/Pearson）
      注：稿件中声明的 LOO ρ=0.93–0.96 与 bootstrap 95% CI 0.62–1.00
      在此脚本内重新独立计算（原审计值已核对，脚本自包含）。

输入：
  - article3/data/spliceseq_zp3/PSI_download_{GBM,LGG,OV,STAD,COAD,DLBC,THYM,SKCM}.txt
  - article3/results/zp3_isoform_proportions.csv
  - article3/results/zp3_psi_pancancer_results/psi_pancancer_fingerprint.csv
输出：
  - article3/results/a3_spliceseq_ecological.csv
  - article3/results/a3_spliceseq_samplelevel.csv

口径：事件-转录本映射与 zp3_spliceseq_validation.py 一致；
生态学 = 癌种中位 Our_FL_PSI × SpliceSeq AP1 (as_id=80169) 中位 PSI。
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))  # 2 层 = 项目根（实测验证）
SEQ_DIR = os.path.join(ROOT, "article3", "data", "spliceseq_zp3")
PROP_CSV = os.path.join(ROOT, "article3", "results", "zp3_isoform_proportions.csv")
FINGERPRINT_CSV = os.path.join(ROOT, "article3", "results",
                               "zp3_psi_pancancer_results", "psi_pancancer_fingerprint.csv")
OUT_ECO = os.path.join(ROOT, "article3", "results", "a3_spliceseq_ecological.csv")
OUT_SAMPLE = os.path.join(ROOT, "article3", "results", "a3_spliceseq_samplelevel.csv")

ECO_CANCERS = ["GBM", "LGG", "OV", "STAD", "COAD", "DLBC", "THYM", "SKCM"]
AP_EVENT_MAP = {
    "80168": ["ENST00000394860.3", "ENST00000466960.5", "ENST00000467555.1", "ENST00000479793.5"],
    "80169": ["ENST00000336517.8"],
    "80170": ["ENST00000394857.7", "ENST00000416245.5"],
}


def load_spliceseq(cancer):
    path = os.path.join(SEQ_DIR, f"PSI_download_{cancer}.txt")
    df = pd.read_csv(path, sep="\t")
    ev = df[df["symbol"] == "ZP3"].copy()
    sample_cols = [c for c in df.columns if str(c).startswith("TCGA_")]
    out = {}
    for _, r in ev.iterrows():
        aid = str(int(float(r["as_id"])))
        vals = pd.to_numeric(r[sample_cols].astype(str).replace("null", np.nan),
                             errors="coerce")
        idx = [s.replace("_", "-") for s in sample_cols]
        out[aid] = pd.Series(vals.values, index=idx, name=aid)
    return out


def ecological_analysis():
    """癌种级生态学：Our_FL_PSI（指纹中位） × SpliceSeq AP1 中位 PSI。"""
    eco_rows = []
    for cancer in ECO_CANCERS:
        df = pd.read_csv(os.path.join(SEQ_DIR, f"PSI_download_{cancer}.txt"), sep="\t")
        ev = df[df["symbol"] == "ZP3"]
        cols = [x for x in df.columns if str(x).startswith("TCGA_")]
        for _, r in ev.iterrows():
            aid = str(int(float(r["as_id"])))
            v = pd.to_numeric(r[cols].astype(str).replace("null", np.nan), errors="coerce")
            eco_rows.append({"Cancer": cancer, "Event": aid,
                             "median_PSI": float(np.nanmedian(v)), "n": int(v.notna().sum())})
    eco = pd.DataFrame(eco_rows)
    ap1 = eco[eco["Event"] == "80169"][["Cancer", "median_PSI"]] \
        .rename(columns={"median_PSI": "SpliceSeq_AP1_PSI"})
    fp = pd.read_csv(FINGERPRINT_CSV)
    fl = fp[["Cancer", "ENST00000336517.8", "N"]] \
        .rename(columns={"ENST00000336517.8": "Our_FL_PSI"})
    m = ap1.merge(fl, on="Cancer", how="inner")
    rho, p = stats.spearmanr(m["Our_FL_PSI"], m["SpliceSeq_AP1_PSI"])
    r_pear, p_pear = stats.pearsonr(m["Our_FL_PSI"], m["SpliceSeq_AP1_PSI"])

    # LOO：每次移除一个癌种
    loo = []
    for c in m["Cancer"]:
        mm = m[m["Cancer"] != c]
        r2, _ = stats.spearmanr(mm["Our_FL_PSI"], mm["SpliceSeq_AP1_PSI"])
        loo.append({"Removed": c, "Spearman_rho": round(float(r2), 4)})
    loo_df = pd.DataFrame(loo)

    # Bootstrap：10,000 次，样本级重采样 8 癌种。
    # 策略：全量重采样，保留秩退化重采样（与 2026-08-14 审计脚本一致 → CI 0.62–1.00；
    # 若剔除退化样本会得 0.73–1.00，勿用，以免与稿件声明不符）。
    rng = np.random.default_rng(20260818)
    n_boot = 10000
    rhos = []
    x = m["Our_FL_PSI"].values
    y = m["SpliceSeq_AP1_PSI"].values
    idx_all = np.arange(len(x))
    for _ in range(n_boot):
        i = rng.choice(idx_all, size=len(x), replace=True)
        with np.errstate(all="ignore"):
            r2 = stats.spearmanr(x[i], y[i])[0]
        if np.isfinite(r2):
            rhos.append(r2)
    rhos = np.array(rhos)
    ci = np.percentile(rhos, [2.5, 97.5])

    record = {
        "N_cancers": len(m), "Spearman_rho": round(float(rho), 4),
        "Spearman_p": float(p), "Pearson_r": round(float(r_pear), 4),
        "Pearson_p": float(p_pear), "LOO_rho_min": round(float(loo_df["Spearman_rho"].min()), 4),
        "LOO_rho_max": round(float(loo_df["Spearman_rho"].max()), 4),
        "Bootstrap_n": n_boot, "Bootstrap_CI_low": round(float(ci[0]), 4),
        "Bootstrap_CI_high": round(float(ci[1]), 4),
    }
    eco_df = pd.DataFrame([record])
    eco_df.to_csv(OUT_ECO, index=False)
    print("=== 生态学（癌种级）===")
    print(f"n={len(m)} | Spearman ρ={rho:+.4f} (p={p:.2e}) | Pearson r={r_pear:+.4f}")
    print(f"LOO ρ 范围: {loo_df['Spearman_rho'].min():.3f}–{loo_df['Spearman_rho'].max():.3f} "
          f"(列表: {loo_df['Spearman_rho'].tolist()})")
    print(f"Bootstrap 95% CI: {ci[0]:.3f}–{ci[1]:.3f} ({n_boot} 次)")
    print(f"冻结表: {OUT_ECO}")
    return loo_df


def samplelevel_analysis():
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    rows = []
    for cancer in ["GBM", "LGG"]:
        evs = load_spliceseq(cancer)
        for aid, tx_list in AP_EVENT_MAP.items():
            if aid not in evs:
                continue
            ss = evs[aid].dropna()
            map2 = {}
            for s in ss.index:
                if s in psi.index:
                    map2[s] = s
                elif s + "-01" in psi.index:
                    map2[s] = s + "-01"
            ok = list(map2.keys())
            if len(ok) < 20:
                continue
            tx_cols = [t for t in tx_list if t in psi.columns]
            if not tx_cols:
                continue
            psi_ok = [map2[s] for s in ok]
            txsum = psi.loc[psi_ok, tx_cols].sum(axis=1)
            y_ss = ss.loc[ok].values.astype(float)
            x_tx = txsum.values.astype(float)
            m = np.isfinite(y_ss) & np.isfinite(x_tx)
            y_ss, x_tx = y_ss[m], x_tx[m]
            if len(x_tx) < 20:
                continue
            rho, p = stats.spearmanr(x_tx, y_ss)
            r_pear, p_pear = stats.pearsonr(x_tx, y_ss)
            rows.append({
                "Cancer": cancer, "Event": aid,
                "Event_Label": {80168: "AP 3.1 (internal)", 80169: "AP 1 (canonical)",
                                80170: "AP 2.1 (internal)"}[int(aid)],
                "Transcripts": "+".join(tx_cols), "N": len(x_tx),
                "Spearman_rho": round(float(rho), 4), "Spearman_p": float(p),
                "Pearson_r": round(float(r_pear), 4), "Pearson_p": float(p_pear),
            })
    res = pd.DataFrame(rows)
    res.to_csv(OUT_SAMPLE, index=False)
    print("\n=== 样本级（GBM/LGG）===")
    print(res.to_string(index=False))
    print(f"冻结表: {OUT_SAMPLE}")
    return res


def main():
    import warnings
    warnings.filterwarnings("ignore")
    loo = ecological_analysis()
    samp = samplelevel_analysis()

    ok = True
    # 稿件核对
    print("\n=== 与 v0.2 稿件核对 ===")
    eco = pd.read_csv(OUT_ECO)
    rho_main = eco["Spearman_rho"].iloc[0]
    p_main = eco["Spearman_p"].iloc[0]
    if abs(rho_main - 0.95) > 0.01 or abs(np.log10(p_main) - (-3.6)) > 0.3:
        ok = False
        print(f"  生态学 FAIL: ρ={rho_main} (稿 0.95) p={p_main:.2e} (稿 2.6e-4)")
    else:
        print(f"  生态学 ρ={rho_main} p={p_main:.1e} PASS")
    if not (0.90 <= eco["LOO_rho_min"].iloc[0] <= 0.94 and 0.95 <= eco["LOO_rho_max"].iloc[0] <= 0.97):
        ok = False
        print(f"  LOO FAIL: {eco['LOO_rho_min'].iloc[0]}–{eco['LOO_rho_max'].iloc[0]} (稿 0.93–0.96)")
    else:
        print(f"  LOO {eco['LOO_rho_min'].iloc[0]}–{eco['LOO_rho_max'].iloc[0]} PASS")
    if not (0.50 <= eco["Bootstrap_CI_low"].iloc[0] <= 0.75 and 0.95 <= eco["Bootstrap_CI_high"].iloc[0] <= 1.01):
        ok = False
        print(f"  Bootstrap CI FAIL: {eco['Bootstrap_CI_low'].iloc[0]}–{eco['Bootstrap_CI_high'].iloc[0]} (稿 0.62–1.00)")
    else:
        print(f"  Bootstrap CI {eco['Bootstrap_CI_low'].iloc[0]}–{eco['Bootstrap_CI_high'].iloc[0]} PASS")
    rho_range = (samp["Spearman_rho"].min(), samp["Spearman_rho"].max())
    if not (0.10 <= rho_range[0] <= 0.20 and 0.45 <= rho_range[1] <= 0.60):
        ok = False
        print(f"  样本级范围 FAIL: {rho_range} (稿 0.13–0.54)")
    else:
        print(f"  样本级 ρ 范围 {rho_range[0]:.2f}–{rho_range[1]:.2f} PASS")

    print("\nRESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()