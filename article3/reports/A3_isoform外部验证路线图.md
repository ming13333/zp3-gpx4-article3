# A3 isoform 级外部验证路线图（2026-08-18）

> 目的：回答「外部 isoform 级验证必须等 dbGaP 吗？能否避开或替代？」
> 结论：**可以避开等待**。已有内部 transportability 补强（冻结⑨）+ 三条替代路径；
> 其中「GEO 重分析」与「RT-qPCR」不依赖权限，可并行推进。

---

## 1. 公开资源重新评估（本轮实测，2026-08-18）

| 数据源 | 层级 | 访问 | 实测结论 |
|---|---|---|---|
| CPTAC-3 "Isoform Expression Quantification" | miRNA isomiR | **open** | **陷阱**：2659 文件与 miRNA EQ 文件数完全相同，文件名均 `.mirnaseq.isoforms.` → 非 mRNA 转录本 isoform（探测脚本 `probe_cptac_isoform_availability.py`，PASS） |
| CPTAC-3 Splice Junction Quantification | junction PSI | controlled | 2968 文件，需 dbGaP 授权 |
| CPTAC-3 Gene Expression Quantification | 基因级 | open | 3042 文件，只能测 gene-level（同外部队列） |
| GLASS-US（GDC） | — | — | 无 RNA-seq（仅 WXS） |
| CGGA mRNAseq_693/325 | 基因级 counts | open 下载 | 无 isoform 级 |
| 公开 GBM long-read（ONT/PacBio） | transcript | 检索未见 | 无现成可下载队列 |

**结论**：公开网格中不存在「可直接下载的 mRNA 转录本 isoform 级肿瘤队列表格」。
但**原始 RNA-seq 重分析**与**实验验证**可以绕开这一限制。

---

## 2. 已完成替代（内部 transportability，冻结⑨，commit 待定）

`freeze_a3_transportability.py` → `a3_transportability_frozen.csv`（1062 行）：

| 检验 | 结果 | 含义 |
|---|---|---|
| LOCO（单癌种留出） | M2 ρ +0.124~+0.144（全量 +0.137） | 非单癌种驱动 |
| **L2CO（两癌种同时留出）** | **M2/Myeloid 各 496 对均 0 方向翻转** | 任意两癌种同时去掉仍稳健 |
| **HELDOUT（预设 seed=42，train 22 / val 10 癌种）** | M2 train +0.144 → val +0.114；Myeloid train +0.123 → val +0.070；**val 方向一致且 CI 不跨 0** | 预设留出癌种上独立复现 |

这是「内部 transportability」，不是外部独立验证；外部 gene-level null 仍如实报告（`a3_external_gbm.csv`）。

---

## 3. 替代路径分层（不依赖 dbGaP 的优先执行）

### 路径 A · GEO 原始 RNA-seq 重分析（推荐，无需权限，成本=工具+算力）

- 原理：GEO 上大量 GBM 队列公开 **FASTQ / BAM**，可下载后自行跑
  Salmon/kallisto 获得转录本定量（ENST 级），再算 FL/(FL+RI) 或 log(FL/RI)。
- 上游参考：Ensembl GRCh38 转录本（ZP3 8 个，FLD=ENST00000336517.8，
  RI=ENST00000466960.5，已在仓库 `zp3_isoform_real_quant.py` 使用）。
- 队列候选：GSE113474（24 例，BAM 公开）、GSE77530（32 例，需核对 FASTQ 供货）、
  CGGA raw reads（需 BIGD 申请，中文通道）。
- 分析流程：
  ```text
  FASTQ/BAM → Salmon (--libType A, txome GRCh38)
    → ENST 级 TPM/counts
    → FL 比例 = FL / ΣZP3, RI 比例, log(FL/RI)
    → 7 免疫评分（同一套签名，z-score 共识）
    → 样本级 Spearman（校正 ZP3 总表达后重算）
  ```
- 成本估算：每 30 样本约 1–3 核·日（Salmon 单样本 20–40 分钟）；无权限障碍。
- **注意**：这是「同一分析管线的外部执行」，是真正的 isoform 级外部验证；
  若结果为 null，仍须如实报告（与 gene-level null 一致叙事）。

### 路径 B · FL/RI 特异性 RT-qPCR（低成本实验，可并行）

- 目标：在不依赖批量测序的条件下独立验证「FL 比例升高与免疫抑制特征共存」。
- 引物设计要点（基于 Ensembl 外显子结构）：
  - **FL 特异**：一对引物跨 FL 独有外显子-外显子剪接连接（如外显子 4→5 边界），
    只扩增剪接完成的 FL 转录本；
  - **RI 特异**：一对引物均位于 retained-intron 内部（intron 保留区域），
    只扩增未剪接的 RI 转录本；
  - **总表达**：引物跨 ZP3 共有外显子（如外显子 2），作为 ZP3 总 RNA 对照；
  - **参考基因**：GAPDH / B2M / HPRT1，ΔΔCt 相对定量；
  - **免疫标志物**：同批样品测 CD68 / CD163 / TREM2 / CD8A，与 FL/RI 比例做相关。
- 样本：10–20 例肿瘤组织 RNA（与已有转录组分析的 cohort 无关的独立样本最理想）。
- ⚠️ **边界声明**：本路线图只给出设计策略，**不提供已验证引物序列**——
  具体引物必须在 wet-lab 用 Primer-BLAST、基因组 BLAT 及 qPCR 熔解曲线验证
  特异性（避免扩增其他 ZP3 转录本/基因组 DNA）。

### 路径 C · dbGaP 申请 CPTAC-3 Splice Junction（长期，非阻塞）

- 申请流程（NIH 官方口径，来源：grants.nih.gov 数据获取页）：
  1. PI 需为机构永久雇员（助理/副教授及以上职称或资深研究者），实验室 PI 角色；
  2. 机构需有 eRA Commons 账户（首次需 business office 注册，含 DUNS）；
  3. PI 在 dbGaP 提交项目申请（研究用途声明 RUS + 非技术摘要）；
  4. 机构 Signing Official 审批并会签；
  5. NIH Data Access Committee（DAC）审核（研究用途须符合该数据集 DUL）；
  6. 批准后有效期 1 年，可续期；需年度进展报告。
- 数据集：CPTAC-3（phs001287）或 CPTAC-2（phs000892），GDC splice junction 受控文件。
- 时间：数周–数月；适合投稿后 revision 阶段补强，不适合作投稿阻塞条件。
- 一旦获批，用我们已定位的 FL/J1–J5、RI/J6–J8 边界算 junction PSI。

### 路径 D · 本地/合作 RNA-seq

- 院内存量 GBM RNA-seq BAM（如有），直接走 Salmon + 本仓库冻结脚本；
- 或合作者提供 quantification 表（ENST 级）。
- 这是最捷径但依赖资源可用性。

### 路径 A 执行状态（2026-08-18 16:40 启动，进行中）

- 队列：**GSE113474 / PRJNA451200**（NYU, Possemato lab; 24 例成人 GBM, 单端 HiSeq2500）。
  与已报告的基因级外部队列 GSE113474 同源——**同一样本, 从 FASTQ 重新做转录本级量化**。
- 工具：**kallisto v0.51.1**（Windows 版, `bin/kallisto/kallisto.exe`）；
  参考转录本 **Ensembl GRCh38 release-110 cdna.all**（1.4 GB）。
- 数据获取：ENA 每样本 ~280 MB FASTQ（24 × ~6.7 GB 总）；网络 ~300 KB/s/连接，
  采用 6 线程 FASTQ 并行 + 8 段 Range 并行下载转录本。
- 分析管线：`article3/scripts/run_external_kallisto_reanalysis.py`（纯标准库 +
  kallisto.exe）：索引 → `kallisto quant --single`（read length 自动检测）→
  ZP3 各转录本 TPM → FL/RI 比例与 log(FL/RI) → 7 免疫评分（z-score 共识，
  同 gene→ENST 映射）→ Spearman → 冻结 `a3_external_isoform_kallisto.csv`。
- 产物目录：`article3/data/external_reanalysis/`（FASTQ/索引 git-ignored，
  脚本与冻结 CSV 入库）。

---

## 4. 决策树

```text
当前投稿版本（v0.3.3）
  需要 isoform 级外部证据吗？
  ├─ 否 → 维持 bounded proxy 叙事（内部多尺度 + 外部 gene-level null + transportability）
  │          ↑ 当前论文采用此路线
  └─ 是 → 走哪条？
        ├─ 有权限制队 BAM/FASTQ → 路径 A（GEO 重分析）→ 本仓库脚本直接算
        ├─ 有实验条件 → 路径 B（RT-qPCR）→ 独立验证 FL/RI 开关
        ├─ 有 PI 资质+机构 eRA → 路径 C（dbGaP）→ 投稿后补强
        └─ 有院内存量测序 → 路径 D → 最快
```

**推荐**：投稿维持「bounded proxy」叙事即可；修订期优先路径 A 或 B。

---

## 5. 对稿件的影响（v0.3.3 已落实）

- Results 新增「Internal cross-cancer transportability」段（LOCO/L2CO/HELDOUT，
  引用 `a3_transportability_frozen.csv`）；
- Limitations 明确「外部 isoform 级验证路径已文档化（GEO 重分析 / RT-qPCR /
  dbGaP junction PSI），等待资源可用」；
- 数据可用性/代码：探测脚本 `probe_cptac_isoform_availability.py` 与路线图文档
  随仓库发布，佐证「已系统排查公开 isoform 数据」。