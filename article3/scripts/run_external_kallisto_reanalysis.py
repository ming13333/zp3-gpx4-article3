#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A3 冻结⑩ — 外部 isoform 级重分析（GEO GSE113474/PRJNA451200, kallisto）
================================================================
目的：在 TCGA/TARGET/GTEx 之外，对 24 例成人 GBM（GSE113474, NYU; Possemato lab,
      PRJNA451200, 单端 HiSeq2500）的公开 FASTQ 做转录本级重分析（kallisto,
      Ensembl GRCh38 cdna release-110），计算 ZP3 FL/RI 比例并检验
      「isoform 比例 × 免疫评分」关联——这是 A3 首个真正的 isoform 级外部验证。

管线（纯标准库 + kallisto.exe）：
  S1. 建索引: kallisto index -i txome/transcripts_k31.idx txome/Homo_sapiens.GRCh38.cdna.all.fa.gz
      （若索引已存在则跳过）
  S2. 量化: kallisto quant -i idx --single -l <read_len> -s 15 -t 4 -o out/<run> fastq/<run>.fastq.gz
      （read_len 从每个样本第一条 read 自动检测；-t 4 并行，4 样本同时跑）
  S3. 解析: 每样本 abundance.tsv → ZP3 各转录本 TPM + 免疫基因 TPM（z-score 共识评分）
  S4. 关联: Spearman(ZP3 FL 比例 / RI 比例 / log(FL/RI), 7 免疫评分)
  S5. 冻结: article3/results/a3_external_isoform_kallisto.csv
  S6. 自检: 阳性 QC（CD8A↔Cytolytic、CD68↔CD163 应为正）+ ZP3 检出 + 方向对照

输入：article3/data/external_reanalysis/fastq/*.fastq.gz（24 样本, 由 download_fastq.sh 获取）
输出：article3/results/a3_external_isoform_kallisto.csv
"""
import os
import sys
import csv
import gzip
import math
import subprocess
import glob

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(BASE))          # 项目根
RE = os.path.join(ROOT, "article3", "data", "external_reanalysis")
FASTQ_DIR = os.path.join(RE, "fastq")
TXOME_GZ = os.path.join(RE, "txome", "Homo_sapiens.GRCh38.cdna.all.fa.gz")
IDX = os.path.join(RE, "txome", "transcripts_k31.idx")
OUT_DIR = os.path.join(RE, "out")
OUT_CSV = os.path.join(ROOT, "article3", "results", "a3_external_isoform_kallisto.csv")
KALLISTO = os.path.join(RE, "bin", "kallisto", "kallisto.exe")

# ZP3 转录本（Ensembl GRCh38, 已在 zp3_isoform_real_quant.py 使用）
ZP3_ISOFORMS = {
    "ENST00000336517.8": "FL",     # 经典全长（肿瘤中比例升高）
    "ENST00000466960.5": "RI",     # retained intron
    "ENST00000394860.3": "T5",     # 5-exon truncated
    "ENST00000467555.1": "T4",
    "ENST00000394857.7": "T3",
    "ENST00000416245.5": "T2",
    "ENST00000479793.5": "T1",
}

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
QC_MARKERS = ["CD8A", "CD68", "CD163"]


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


# ---------------------------------------------------------------------------
def read_len_of(fastq_gz):
    """读第一条 read 的序列长度（单端）"""
    with gzip.open(fastq_gz, "rt") as f:
        f.readline()  # @header
        seq = f.readline().strip()
        return len(seq)


def build_index():
    if os.path.isfile(IDX):
        print(f"索引已存在: {IDX}"); return
    assert os.path.isfile(TXOME_GZ), f"缺失转录本: {TXOME_GZ}"
    print("建索引 (k31)...")
    subprocess.run([KALLISTO, "index", "-i", IDX, TXOME_GZ], check=True)
    print("索引完成")


def parse_abundance(path):
    """返回 {ENST(去版本): tpm}"""
    out = {}
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for d in rd:
            tid = d["target_id"].split(".")[0]
            out[tid] = f2(d["tpm"])
    return out


def detect_fastqs():
    fs = sorted(glob.glob(os.path.join(FASTQ_DIR, "*.fastq.gz")))
    return [f for f in fs if os.path.getsize(f) > 1e6]


def check_inputs():
    """下载就绪检查：FASTQ ≥24 且完整、txome 完整；不满足则提示进度并返回 False。"""
    missing = []
    fastqs = detect_fastqs()
    n_ok = 0
    for f in fastqs:
        try:
            with gzip.open(f, "rb") as g:
                while g.read(1 << 20):
                    pass
            n_ok += 1
        except Exception:
            missing.append(os.path.basename(f))
    tx_ok = os.path.isfile(TXOME_GZ) and os.path.getsize(TXOME_GZ) > 800_000_000
    print(f"[check] FASTQ 完整 {n_ok}/24 (不完整 {len(missing)}), txome {'OK' if tx_ok else '未就绪/不完整'}")
    if not tx_ok:
        print("  txome 未完成: 等 download_txome.py 结束")
    if n_ok < 24:
        print(f"  FASTQ 未齐: 等 download_fastq.py 结束 (当前 {n_ok}/24)")
    return n_ok >= 24 and tx_ok


def main():
    assert os.path.isfile(KALLISTO), f"缺失 kallisto: {KALLISTO}"
    if not check_inputs():
        print("下载未完成，请稍后重跑本脚本（自动进入索引→量化→分析）")
        return 2
    fastqs = detect_fastqs()
    print(f"检测到 FASTQ 样本: {len(fastqs)}")

    build_index()

    # S2: 量化（4 路并行）
    os.makedirs(OUT_DIR, exist_ok=True)
    todo = []
    for fq in fastqs:
        run = os.path.basename(fq).replace(".fastq.gz", "")
        out_run = os.path.join(OUT_DIR, run, "abundance.tsv")
        if not os.path.isfile(out_run):
            todo.append((run, fq))
    print(f"待量化: {len(todo)} / 总数 {len(fastqs)}")
    procs = []
    for run, fq in todo:
        rl = read_len_of(fq)
        cmd = [KALLISTO, "quant", "-i", IDX, "--single", "-l", str(rl),
               "-s", "15", "-t", "3", "-o", os.path.join(OUT_DIR, run), fq]
        print(f"  quant {run} (read_len={rl})")
        p = subprocess.Popen(cmd)
        procs.append((run, p))
        while len([x for _, x in procs if x.poll() is None]) >= 4:
            import time
            time.sleep(5)
    for run, p in procs:
        p.wait()
        if p.returncode != 0:
            print(f"  FAIL quant {run}")
            return 1
    print("量化完成")

    # S3: 汇总样本级矩阵
    samples = []
    for fq in fastqs:
        run = os.path.basename(fq).replace(".fastq.gz", "")
        ab = parse_abundance(os.path.join(OUT_DIR, run, "abundance.tsv"))
        samples.append((run, ab))
    n = len(samples)
    print(f"样本数: {n}")

    # ZP3 isoform TPM 矩阵
    zp3_rows = {}
    for tid in ZP3_ISOFORMS:
        zp3_rows[tid] = [ab.get(tid, 0.0) for _, ab in samples]
    zp3_total = [sum(ab.get(t, 0.0) for t in ZP3_ISOFORMS) for _, ab in samples]

    def prop(tid):
        return [ (zp3_rows[tid][i] / zp3_total[i]) if zp3_total[i] > 0 else float("nan")
                 for i in range(n) ]
    fl = prop("ENST00000336517")
    ri = prop("ENST00000466960")
    log_ratio = [math.log((fl[i] + 1e-9) / (ri[i] + 1e-9)) for i in range(n)]
    print(f"ZP3 检出: {sum(1 for t in zp3_total if t > 0)}/{n} 样本; "
          f"FL 比例范围 [{min(fl):.3f}, {max(fl):.3f}]")

    # gene symbol → ENST 映射（从 cdna header 解析，缓存）
    sym_map = {}
    cache = os.path.join(RE, "txome", "sym2enst.json")
    if os.path.isfile(cache):
        import json
        sym_map = json.load(open(cache))
    else:
        print("解析 cdna header 建 gene→ENST 映射...")
        with gzip.open(TXOME_GZ, "rt") as f:
            cur = None
            for line in f:
                if line.startswith(">"):
                    h = line[1:].strip().split()
                    tid = h[0].split(".")[0]
                    sym = None
                    for part in h[1:]:
                        if part.startswith("gene:"):
                            sym = part.split(":", 1)[1]
                    if sym:
                        sym_map.setdefault(sym, []).append(tid)
                    cur = sym
        import json
        json.dump(sym_map, open(cache, "w"))
        print(f"映射基因数: {len(sym_map)}")

    # 免疫评分（z-score 共识）
    def gene_score(symbols):
        tids = []
        for s in symbols:
            tids.extend(sym_map.get(s, []))
        # 聚合到基因 TPM
        gene_tpm = {}
        for tid in tids:
            for i, (run, ab) in enumerate(samples):
                v = ab.get(tid, 0.0)
                gene_tpm.setdefault(i, 0.0)
                gene_tpm[i] += v  # 转录本 TPM 求和近似基因 TPM
        vals = [gene_tpm[i] for i in range(n)]
        m = sum(vals) / n
        sd = math.sqrt(sum((v - m) ** 2 for v in vals) / n) or 1.0
        return [(v - m) / sd for v in vals]

    # 基因级 TPM（免疫基因 + QC 标记）统一取各基因 TPM
    scores = {}
    for feat, genes in IMMUNE.items():
        scores[feat] = gene_score(genes)
    # QC 标记单基因 TPM（CD8A/CD68/CD163 及 cytolytic 基因）
    def gene_tpm(symbol):
        tids = sym_map.get(symbol, [])
        out = []
        for run, ab in samples:
            out.append(sum(ab.get(t, 0.0) for t in tids))
        return out

    cd8a = gene_tpm("CD8A")
    cd68 = gene_tpm("CD68")
    cd163 = gene_tpm("CD163")
    cyt = gene_score(["GZMA", "GZMB", "PRF1", "IFNG"])

    # S4: Spearman 关联
    rows_out = []
    print("\n=== ZP3 isoform 比例 × 免疫评分 (外部 GBM, kallisto) ===")
    for feat in IMMUNE:
        r_fl, p_fl = spearman(fl, scores[feat])
        r_ri, p_ri = spearman(ri, scores[feat])
        r_lr, p_lr = spearman(log_ratio, scores[feat])
        print(f"  {feat:18s} FL ρ={r_fl:+.3f} (p={p_fl:.2e}) | "
              f"RI ρ={r_ri:+.3f} | log(FL/RI) ρ={r_lr:+.3f}")
        rows_out.append({
            "Cohort": "GSE113474_kallisto", "N": n, "Metric": "FL_proportion",
            "Feature": feat, "rho": round(r_fl, 4), "p": f"{p_fl:.3e}",
            "Direction_vs_TCGA": "same" if r_fl > 0 else "opposite",
        })
        rows_out.append({
            "Cohort": "GSE113474_kallisto", "N": n, "Metric": "RI_proportion",
            "Feature": feat, "rho": round(r_ri, 4), "p": f"{p_ri:.3e}",
            "Direction_vs_TCGA": "same" if r_ri < 0 else "opposite",
        })
        rows_out.append({
            "Cohort": "GSE113474_kallisto", "N": n, "Metric": "log_FL_over_RI",
            "Feature": feat, "rho": round(r_lr, 4), "p": f"{p_lr:.3e}",
            "Direction_vs_TCGA": "same" if r_lr > 0 else "opposite",
        })

    # S5: 冻结
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"\n冻结表: {OUT_CSV}")

    # S6: 自检
    print("\n=== 自检 ===")
    ok = True
    r1, _ = spearman(cd8a, cyt)
    r2, _ = spearman(cd68, cd163)
    print(f"  QC Spearman(CD8A, Cytolytic)={r1:+.3f}  Spearman(CD68, CD163)={r2:+.3f}")
    if not (r1 > 0 and r2 > 0):
        print("  FAIL: 免疫评分 QC 未通过"); ok = False
    else:
        print("  PASS: 免疫评分有效")
    fl_nonzero = sum(1 for t in zp3_total if t > 0)
    if fl_nonzero < n * 0.8:
        print(f"  WARN: ZP3 检出率低 ({fl_nonzero}/{n})")
    # 方向对照 TCGA 内部：M2 FL 正、Myeloid FL 正
    for feat in ["M2_Macrophage", "Myeloid"]:
        row = [r for r in rows_out if r["Metric"] == "FL_proportion" and r["Feature"] == feat][0]
        print(f"  {feat}: 外部 FL ρ={row['rho']:+.3f} ({row['Direction_vs_TCGA']} vs TCGA 内部 +)")
        if row["Direction_vs_TCGA"] != "same":
            print(f"  WARN: {feat} 方向与 TCGA 内部不一致（如实报告，不判失败）")
    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def IMMUNE_genes():
    return sorted({g for v in IMMUNE.values() for g in v})


if __name__ == "__main__":
    sys.exit(main())