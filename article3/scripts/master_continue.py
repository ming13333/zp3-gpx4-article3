#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
master_continue.py — Path A Stage ② master control (run once after FASTQ are ready)
================================================================
Prerequisite: the 24 FASTQ files under article3/data/external_reanalysis/fastq/ are complete (download_fastq.py done)
Functions:
  1. Check 24 FASTQ completeness + txome completeness; if either is missing, fill it first
  2. Transcriptome download (download_txome.py, 8-segment Range + retry) → gzip verification
  3. Call run_external_kallisto_reanalysis.py (index → 24-sample quantification → ZP3 FL/RI →
     7 immune scores → Spearman → freeze a3_external_isoform_kallisto.csv → self-check)
Usage: python article3/scripts/master_continue.py
"""
import os
import sys
import gzip
import glob
import subprocess

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))
RE = os.path.join(ROOT, "article3", "data", "external_reanalysis")
FASTQ_DIR = os.path.join(RE, "fastq")
TXOME_GZ = os.path.join(RE, "txome", "Homo_sapiens.GRCh38.cdna.all.fa.gz")
PY = sys.executable


def check_fastqs():
    fs = sorted(glob.glob(os.path.join(FASTQ_DIR, "*.fastq.gz")))
    ok = []
    for f in fs:
        try:
            with gzip.open(f, "rb") as g:
                while g.read(1 << 20):
                    pass
            ok.append(os.path.basename(f))
        except Exception:
            pass
    return ok


def check_txome():
    if not os.path.isfile(TXOME_GZ):
        return False
    try:
        with gzip.open(TXOME_GZ, "rb") as g:
            while g.read(1 << 20):
                pass
        return os.path.getsize(TXOME_GZ) > 800_000_000
    except Exception:
        return False


def main():
    fq_ok = check_fastqs()
    print(f"[1] FASTQ complete {len(fq_ok)}/24")
    if len(fq_ok) < 24:
        print("    FASTQ not all present (download_fastq.py not yet complete); this run only reports progress.")
        print("    Please wait until FASTQ download is complete and then rerun this script.")
        return 0 if len(fq_ok) >= 8 else 2

    tx_ok = check_txome()
    print(f"[2] txome complete: {tx_ok}")
    if not tx_ok:
        print("    Starting download_txome.py (transcriptome 8-segment download)...")
        r = subprocess.run([PY, os.path.join(RE, "download_txome.py")])
        if r.returncode != 0 or not check_txome():
            print("FAIL: transcriptome download/verification failed")
            return 1
        print("    transcriptome ready")

    print("[3] Running full reanalysis pipeline (index→quantification→analysis→freeze)...")
    r = subprocess.run([PY, os.path.join(BASE, "run_external_kallisto_reanalysis.py")])
    print(f"[3] pipeline exit code: {r.returncode}")
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
