# -*- coding: utf-8 -*-
"""
Article 3 Supplement ① — SpliceSeq event-level PSI cross-validation
================================================
Purpose: verify whether "transcript TPM proportion proxy PSI" is consistent with independent event-level PSI (TCGA SpliceSeq，
based on junction/exon read counts, independent of TPM quantification）.

Data sources:
  - spliceseq_zp3/PSI_download_GBM.txt / PSI_download_LGG.txt
    （downloaded via POST https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIDownload，
     ZP3 has 3 AP events: 80168 exons=3.1 / 80169 exons=1 / 80170 exons=2.1）
  - zp3_isoform_proportions.csv（19131 samples × 7 transcript proportions）

Event↔transcript mapping（based on Ensembl GRCh38 5' end structure）:
  - AP 80169 (exons=1, most 5' classical promoter)  ↔ FL canonical ENST00000336517.8
  - AP 80170 (exons=2.1, internal promoter)      ↔ ENST00000394857/00416245（intermediate transcript）
  - AP 80168 (exons=3.1, most 3' internal promoter) ↔ ENST00000394860/00466960(RI)/00467555/00479793

Criterion: Spearman ρ > 0.7 → proxy PSI methodology stands; otherwise downgraded to TRA definition.
Output: spliceseq_zp3/spliceseq_validation_results.csv + fig_spliceseq_validation.png
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

# Event → transcript mapping
AP_EVENT_MAP = {
    "80168": {  # exons=3.1, most 3' internal promoter
        "label": "AP 3.1 (internal promoter, 3')",
        "txs": ["ENST00000394860.3", "ENST00000466960.5", "ENST00000467555.1", "ENST00000479793.5"],
    },
    "80169": {  # exons=1, canonical promoter
        "label": "AP 1 (canonical promoter, 5')",
        "txs": ["ENST00000336517.8"],  # FL canonical
    },
    "80170": {  # exons=2.1, internal promoter
        "label": "AP 2.1 (internal promoter)",
        "txs": ["ENST00000394857.7", "ENST00000416245.5"],
    },
}

# Cancers for ecological validation (FL fingerprint high/low extremes + glioma)
ECO_CANCERS = ["GBM", "LGG", "OV", "STAD", "COAD", "DLBC", "THYM", "SKCM"]


def load_spliceseq(cancer):
    """Read SpliceSeq PSI file, return event_id -> (label, Series[sample->PSI])."""
    path = os.path.join(SEQ_DIR, f"PSI_download_{cancer}.txt")
    df = pd.read_csv(path, sep="\t")
    ev = df[df["symbol"] == "ZP3"].copy()
    sample_cols = [c for c in df.columns if str(c).startswith("TCGA_")]
    out = {}
    for _, r in ev.iterrows():
        aid = str(int(float(r["as_id"])))
        vals = r[sample_cols].astype(str).replace("null", np.nan).replace("", np.nan)
        vals = pd.to_numeric(vals, errors="coerce")
        # SpliceSeq sample name TCGA_02_0047 -> TCGA-02-0047 (aligned with proportion matrix prefix)
        idx = [s.replace("_", "-") for s in sample_cols]
        s = pd.Series(vals.values, index=idx, name=aid)
        out[aid] = s
    return out


def load_proportions():
    psi = pd.read_csv(PROP_CSV, index_col=0)
    psi.columns = [c.strip() for c in psi.columns]
    return psi


def main():
    print("=== Article 3 Supplement ①: SpliceSeq event-level PSI cross-validation ===\n")
    psi = load_proportions()
    print(f"Proportion matrix: {psi.shape[0]} samples × {psi.shape[1]} transcripts")

    rows = []
    panels = []
    for cancer in ["GBM", "LGG"]:
        evs = load_spliceseq(cancer)
        print(f"\n[{cancer}] SpliceSeq events: {list(evs.keys())}")
        for aid, meta in AP_EVENT_MAP.items():
            if aid not in evs:
                print(f"  !! Event {aid} missing")
                continue
            ss = evs[aid].dropna()
            # Sample prefix alignment: SpliceSeq uses TCGA-XX-XXXX, ratio matrix uses TCGA-XX-XXXX-01
            # Build ss sample -> psi sample mapping
            map2 = {}
            for s in ss.index:
                if s in psi.index:
                    map2[s] = s
                elif s + "-01" in psi.index:
                    map2[s] = s + "-01"
            ok = list(map2.keys())
            if len(ok) < 20:
                print(f"  {meta['label']} (as_id={aid}): only {len(ok)} samples aligned, skipping")
                continue
            # Aggregate transcript proportions for event (internal promoter = sum of multiple transcripts)
            tx_cols = [t for t in meta["txs"] if t in psi.columns]
            if not tx_cols:
                print(f"  {meta['label']}: no corresponding transcript columns")
                continue
            psi_ok = [map2[s] for s in ok]
            txsum = psi.loc[psi_ok, tx_cols].sum(axis=1)
            y_ss = ss.loc[ok].values.astype(float)
            x_tx = txsum.values.astype(float)
            m = np.isfinite(y_ss) & np.isfinite(x_tx)
            y_ss = y_ss[m]
            x_tx = x_tx[m]
            if len(x_tx) < 20:
                print(f"  {meta['label']}: valid samples {len(x_tx)} < 20, skipping")
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
    print(f"\nResults saved: {OUT_CSV}")

    # ---- Figure ----
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
        print(f"Figure saved: {OUT_FIG}")

    # ---- Cross-cancer ecological validation (sample-level → cancer-level) ----
    print("\n=== Cross-cancer ecological validation ===")
    eco_rows = []
    for cancer in ECO_CANCERS:
        path = os.path.join(SEQ_DIR, f"PSI_download_{cancer}.txt")
        if not os.path.exists(path):
            print(f"  !! {cancer} data missing, skipping")
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
        print(f"  Cancer-level n={len(m)}: Our FL-PSI × SpliceSeq AP-1 PSI")
        print(f"    Spearman ρ={rho_e:+.3f} (p={p_e:.2e}) | Pearson r={rho_pe:+.3f} (p={p_pe:.2e})")
        m.to_csv(os.path.join(SEQ_DIR, "spliceseq_ecological_validation.csv"), index=False)
        # Ecological scatter plot
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
        print(f"  Ecological plot: fig_spliceseq_ecological.png")
    else:
        print(f"  Insufficient cancer matching ({len(m)})")

    # ---- Verdict ----
    print("\n=== Verdict ===")
    if len(res) == 0:
        print("!! No valid validation results")
    for _, r in res.iterrows():
        verdict = "✅ Holds (ρ>0.7)" if r["Spearman_rho"] > 0.7 else (
            "⚠️ Moderate (0.4<ρ≤0.7)" if r["Spearman_rho"] > 0.4 else "❌ Inconsistent")
        print(f"  {r['Cancer']} {r['Event_Label']}: ρ={r['Spearman_rho']:+.3f} → {verdict}")
    if len(m) >= 5:
        v2 = "✅ Holds (ρ>0.7)" if rho_e > 0.7 else "⚠️ Threshold not met"
        print(f"  Ecological (cancer-level): ρ={rho_e:+.3f} → {v2}")


if __name__ == "__main__":
    main()
