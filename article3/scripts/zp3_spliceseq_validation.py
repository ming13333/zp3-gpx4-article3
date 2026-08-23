# -*- coding: utf-8 -*-
"""
Article 3 补强① — SpliceSeq 事件级 PSI 交叉验证
================================================
目的：验证"转录本 TPM 比例代理 PSI"是否与独立的事件级 PSI（TCGA SpliceSeq，
基于 junction/exon 读取计数，与 TPM 定量相互独立）一致。

数据源：
  - spliceseq_zp3/PSI_download_GBM.txt / PSI_download_LGG.txt
    （POST https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIDownload 下载，
     ZP3 共 3 个 AP 事件：80168 exons=3.1 / 80169 exons=1 / 80170 exons=2.1）
  - zp3_isoform_proportions.csv（19131 样本 × 7 转录本比例）

事件↔转录本映射（依据 Ensembl GRCh38 5' 端结构）：
  - AP 80169 (exons=1, 最 5' 经典启动子)  ↔ FL canonical ENST00000336517.8
  - AP 80170 (exons=2.1, 内部启动子)      ↔ ENST00000394857/00416245（中间转录本）
  - AP 80168 (exons=3.1, 最 3' 内部启动子) ↔ ENST00000394860/00466960(RI)/00467555/00479793

判定：Spearman ρ > 0.7 → 代理 PSI 方法学站住；否则降级为 TRA 口径。
产物：spliceseq_zp3/spliceseq_validation_results.csv + fig_spliceseq_validation.png
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "data", "spliceseq_zp3")
PROP_CSV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_isoform_proportions.csv")
OUT_CSV = os.path.join(SEQ_DIR, "spliceseq_validation_results.csv")
OUT_FIG = os.path.join(SEQ_DIR, "fig_spliceseq_validation.png")

# 事件 → 转录本 映射
AP_EVENT_MAP = {
    "80168": {  # exons=3.1, 最 3' 内部启动子
        "label": "AP 3.1 (internal promoter, 3')",
        "txs": ["ENST00000394860.3", "ENST00000466960.5", "ENST00000467555.1", "ENST00000479793.5"],
    },
    "80169": {  # exons=1, 经典启动子
        "label": "AP 1 (canonical promoter, 5')",
        "txs": ["ENST00000336517.8"],  # FL canonical
    },
    "80170": {  # exons=2.1, 中间启动子
        "label": "AP 2.1 (internal promoter)",
        "txs": ["ENST00000394857.7", "ENST00000416245.5"],
    },
}

# 生态学验证的癌种（FL 指纹高/低两端 + 胶质瘤）
ECO_CANCERS = ["GBM", "LGG", "OV", "STAD", "COAD", "DLBC", "THYM", "SKCM"]


def load_spliceseq(cancer):
    """读取 SpliceSeq PSI 文件，返回 event_id -> (label, Series[样本->PSI])。"""
    path = os.path.join(SEQ_DIR, f"PSI_download_{cancer}.txt")
    df = pd.read_csv(path, sep="\t")
    ev = df[df["symbol"] == "ZP3"].copy()
    sample_cols = [c for c in df.columns if str(c).startswith("TCGA_")]
    out = {}
    for _, r in ev.iterrows():
        aid = str(int(float(r["as_id"])))
        vals = r[sample_cols].astype(str).replace("null", np.nan).replace("", np.nan)
        vals = pd.to_numeric(vals, errors="coerce")
        # SpliceSeq 样本名 TCGA_02_0047 -> TCGA-02-0047（与比例矩阵前缀对齐）
        idx = [s.replace("_", "-") for s in sample_cols]
        s = pd.Series(vals.values, index=idx, name=aid)
        out[aid] = s
    return out


def load_proportions():
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    return psi


def main():
    print("=== Article 3 补强①: SpliceSeq 事件级 PSI 交叉验证 ===\n")
    psi = load_proportions()
    print(f"比例矩阵: {psi.shape[0]} 样本 × {psi.shape[1]} 转录本")

    rows = []
    panels = []
    for cancer in ["GBM", "LGG"]:
        evs = load_spliceseq(cancer)
        print(f"\n[{cancer}] SpliceSeq 事件: {list(evs.keys())}")
        for aid, meta in AP_EVENT_MAP.items():
            if aid not in evs:
                print(f"  !! 事件 {aid} 缺失")
                continue
            ss = evs[aid].dropna()
            # 样本前缀对齐：SpliceSeq 用 TCGA-XX-XXXX，比例矩阵用 TCGA-XX-XXXX-01
            # 构建 ss 样本 -> psi 样本 映射
            map2 = {}
            for s in ss.index:
                if s in psi.index:
                    map2[s] = s
                elif s + "-01" in psi.index:
                    map2[s] = s + "-01"
            ok = list(map2.keys())
            if len(ok) < 20:
                print(f"  {meta['label']} (as_id={aid}): 仅 {len(ok)} 样本对齐，跳过")
                continue
            # 聚合事件对应的转录本比例（内部启动子 = 多个转录本之和）
            tx_cols = [t for t in meta["txs"] if t in psi.columns]
            if not tx_cols:
                print(f"  {meta['label']}: 无对应转录本列")
                continue
            psi_ok = [map2[s] for s in ok]
            txsum = psi.loc[psi_ok, tx_cols].sum(axis=1)
            y_ss = ss.loc[ok].values.astype(float)
            x_tx = txsum.values.astype(float)
            m = np.isfinite(y_ss) & np.isfinite(x_tx)
            y_ss = y_ss[m]
            x_tx = x_tx[m]
            if len(x_tx) < 20:
                print(f"  {meta['label']}: 有效样本 {len(x_tx)} < 20，跳过")
                continue
            rho, p = stats.spearmanr(x_tx, y_ss)
            rho_pear, p_pear = stats.pearsonr(x_tx, y_ss)
            rows.append({
                "Cancer": cancer, "Event": aid, "Event_Label": meta["label"],
                "Transcripts": "+".join(tx_cols), "N": len(x_tx),
                "Spearman_rho": round(float(rho), 4), "Spearman_p": float(p),
                "Pearson_r": round(float(rho_pear), 4), "Pearson_p": float(p_pear),
            })
            panels.append((cancer, meta["label"], aid, x_tx, y_ss, rho, p, len(x_tx)))
            print(f"  {meta['label']} (as_id={aid}): n={len(x_tx)} Spearman ρ={rho:+.3f} (p={p:.2e}) Pearson r={rho_pear:+.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(OUT_CSV, index=False)
    print(f"\n结果已存: {OUT_CSV}")

    # ---- 图 ----
    if panels:
        n = len(panels)
        fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 4.6))
        if n == 1:
            axes = [axes]
        for ax, (cancer, label, aid, x, y, rho, p, nn) in zip(axes, panels):
            ax.scatter(x, y, s=30, alpha=0.6, c="#378ADD", edgecolor="white", linewidth=0.4)
            z = np.polyfit(x, y, 1)
            xs = np.linspace(np.nanmin(x), np.nanmax(x), 50)
            ax.plot(xs, np.polyval(z, xs), "r--", lw=1.2, alpha=0.7)
            ax.set_xlabel("Our transcript proportion (TPM-ratio)", fontsize=10)
            ax.set_ylabel("SpliceSeq event PSI", fontsize=10)
            ax.set_title(f"{cancer} | {label}\nn={nn}, Spearman ρ={rho:+.3f}, p={p:.2e}",
                         fontsize=10)
            ax.axvline(0, color="grey", lw=0.5, ls=":")
            ax.axhline(0, color="grey", lw=0.5, ls=":")
        fig.suptitle("Cross-validation: TPM-ratio proxy vs SpliceSeq event-level PSI (ZP3)",
                     fontsize=13, y=1.02)
        fig.tight_layout()
        fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"图已存: {OUT_FIG}")

    # ---- 跨癌种生态学验证（样本级 → 癌种级）----
    print("\n=== 跨癌种生态学验证 ===")
    eco_rows = []
    for cancer in ECO_CANCERS:
        path = os.path.join(SEQ_DIR, f"PSI_download_{cancer}.txt")
        if not os.path.exists(path):
            print(f"  !! {cancer} 数据缺失，跳过")
            continue
        df = pd.read_csv(path, sep="\t")
        ev = df[df["symbol"] == "ZP3"]
        cols = [x for x in df.columns if str(x).startswith("TCGA_")]
        for _, r in ev.iterrows():
            aid = str(int(float(r["as_id"])))
            v = pd.to_numeric(r[cols].astype(str).replace("null", np.nan), errors="coerce")
            eco_rows.append({"Cancer": cancer, "Event": aid,
                             "median_PSI": float(np.nanmedian(v)), "n": int(v.notna().sum())})
    eco = pd.DataFrame(eco_rows)
    eco.to_csv(os.path.join(SEQ_DIR, "spliceseq_eco_by_cancer.csv"), index=False)
    ap1 = eco[eco["Event"] == "80169"][["Cancer", "median_PSI"]] \
        .rename(columns={"median_PSI": "SpliceSeq_AP1_PSI"})
    fp = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(BASE))), "article3", "results", "zp3_psi_pancancer_results", "psi_pancancer_fingerprint.csv"))
    fl = fp[["Cancer", "ENST00000336517.8", "N"]].rename(columns={"ENST00000336517.8": "Our_FL_PSI"})
    m = ap1.merge(fl, on="Cancer", how="inner")
    if len(m) >= 5:
        rho_e, p_e = stats.spearmanr(m["Our_FL_PSI"], m["SpliceSeq_AP1_PSI"])
        rho_pe, p_pe = stats.pearsonr(m["Our_FL_PSI"], m["SpliceSeq_AP1_PSI"])
        print(f"  癌种级 n={len(m)}: Our FL-PSI × SpliceSeq AP-1 PSI")
        print(f"    Spearman ρ={rho_e:+.3f} (p={p_e:.2e}) | Pearson r={rho_pe:+.3f} (p={p_pe:.2e})")
        m.to_csv(os.path.join(SEQ_DIR, "spliceseq_ecological_validation.csv"), index=False)
        # 生态学散点图
        fig2, ax2 = plt.subplots(figsize=(6.2, 5))
        ax2.scatter(m["Our_FL_PSI"], m["SpliceSeq_AP1_PSI"], s=90, alpha=0.85,
                    c="#C00000", edgecolor="white", zorder=3)
        for _, r in m.iterrows():
            ax2.annotate(r["Cancer"], (r["Our_FL_PSI"], r["SpliceSeq_AP1_PSI"]),
                         textcoords="offset points", xytext=(6, 5), fontsize=9)
        z = np.polyfit(m["Our_FL_PSI"], m["SpliceSeq_AP1_PSI"], 1)
        xs = np.linspace(m["Our_FL_PSI"].min(), m["Our_FL_PSI"].max(), 50)
        ax2.plot(xs, np.polyval(z, xs), "k--", lw=1.2, alpha=0.7)
        ax2.set_xlabel("Our FL-canonical transcript proportion (median per cancer)", fontsize=10)
        ax2.set_ylabel("SpliceSeq AP-1 event PSI (median per cancer)", fontsize=10)
        ax2.set_title(f"Ecological cross-validation (n={len(m)} cancers)\n"
                      f"Spearman ρ={rho_e:+.3f}, p={p_e:.2e}", fontsize=11)
        fig2.tight_layout()
        fig2.savefig(os.path.join(SEQ_DIR, "fig_spliceseq_ecological.png"), dpi=200)
        plt.close(fig2)
        print(f"  生态学图: fig_spliceseq_ecological.png")
    else:
        print(f"  癌种匹配不足 ({len(m)})")

    # ---- 判定 ----
    print("\n=== 判定 ===")
    if len(res) == 0:
        print("!! 无有效验证结果")
    for _, r in res.iterrows():
        verdict = "✅ 站住 (ρ>0.7)" if r["Spearman_rho"] > 0.7 else (
            "⚠️ 中等 (0.4<ρ≤0.7)" if r["Spearman_rho"] > 0.4 else "❌ 不一致")
        print(f"  {r['Cancer']} {r['Event_Label']}: ρ={r['Spearman_rho']:+.3f} → {verdict}")
    if len(m) >= 5:
        v2 = "✅ 站住 (ρ>0.7)" if rho_e > 0.7 else "⚠️ 未达阈值"
        print(f"  生态学(癌种级): ρ={rho_e:+.3f} → {v2}")


if __name__ == "__main__":
    main()
