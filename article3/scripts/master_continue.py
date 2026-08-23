#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
master_continue.py — 路径 A 阶段② 总控（FASTQ 齐备后运行一次）
================================================================
前提：article3/data/external_reanalysis/fastq/ 下 24 个 FASTQ 已完整（download_fastq.py 完成）
功能：
  1. 检查 24 FASTQ 完整 + txome 完整；缺哪个先补哪个
  2. 转录本下载（download_txome.py，8 段 Range + 重试）→ gzip 校验
  3. 调用 run_external_kallisto_reanalysis.py（索引→24 样本量化→ZP3 FL/RI→
     7 免疫评分→Spearman→冻结 a3_external_isoform_kallisto.csv→自检）
用法：python article3/scripts/master_continue.py
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
    print(f"[1] FASTQ 完整 {len(fq_ok)}/24")
    if len(fq_ok) < 24:
        print("    FASTQ 未齐（download_fastq.py 尚未完成），本次仅报告进度。")
        print("    请等待 FASTQ 下载完成后重跑本脚本。")
        return 0 if len(fq_ok) >= 8 else 2

    tx_ok = check_txome()
    print(f"[2] txome 完整: {tx_ok}")
    if not tx_ok:
        print("    启动 download_txome.py（转录本 8 段下载）...")
        r = subprocess.run([PY, os.path.join(RE, "download_txome.py")])
        if r.returncode != 0 or not check_txome():
            print("FAIL: 转录本下载/校验失败")
            return 1
        print("    转录本就绪")

    print("[3] 运行完整重分析管线（索引→量化→分析→冻结）...")
    r = subprocess.run([PY, os.path.join(BASE, "run_external_kallisto_reanalysis.py")])
    print(f"[3] 管线退出码: {r.returncode}")
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())