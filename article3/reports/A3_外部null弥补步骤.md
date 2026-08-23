# A3 外部 null 弥补 · 可执行步骤清单

> 目标：把「外部基因级复制 null」转化为论文的方法学论据，而不是把它救回阳性。
> 原则：不选择性报告、不换单侧检验、不挑队列。所有新分析均冻结进仓库。

---

## 阶段 0 · 先建立记账清单（不动分析）

| # | 动作 | 产出 | 状态 |
|---|---|---|---|
| 0.1 | 在 `article3/results/` 增加 `a3_external_null_diagnostics.md`，记录 null 的约束解释 | 记账文档 | 建议先建 |
| 0.2 | 检查 git 状态，确认冻结表已在版本库 | `git status` | 建议先建 |

> 目的：让「null 是结果而非失败」先被文档化，后续所有分析都在这个框架下展开。

---

## 阶段 1 · 外部队列 null 诊断（最优先，当前数据即可完成）

### 1.1 ZP3 检测质量与分布诊断

- 统计两队列的：
  - ZP3 零值比例；非零样本数；IQR；中位数；检测范围。
- 若 GSE77530 中 ZP3 零值/近零占比高，说明低表达平台无法稳定测量该基因。
- 若大量并列秩（tied ranks），秩相关本身会被压缩。

### 1.2 功效分析（关键数字）

- 对 n=32 和 n=24，计算双侧 Spearman 在 α=0.05 下可检测的最小 |ρ|。
- 输出：「当前队列最多只能检测 |ρ|≥0.35 / 0.40 的关联；0.2–0.3 的真实效应在此样本量下大概率落到 null 区间」。
- 该数字直接写进 Limitations。

### 1.3 表达尺度敏感性

- 对 ZP3 分别做：
  - 原始表达；
  - log1p 变换；
  - 检测/未检测二值；
  - 高表达（≥上四分位）与低表达分层。
- 观察 ρ 方向是否在尺度间稳定。尺度改变后仍为 null，说明不是简单变换问题；
  若尺度间方向翻转，则是测量伪影。

### 1.4 协变量敏感性（若元数据允许）

- 若两队列提供纯度/IDH/年龄/性别/治疗史，加入分层或偏相关校正。
- 若拿不到，如实说明「无法校正 confounder，故 null 仍混杂未校正」。

### 1.5 稳定化检验

- 留一样本法（jackknife）或 bootstrap 重抽样，确认 ρ 不依赖少数样本。
- 排除「单个高表达样本驱动微弱负相关」或「单个样本把正相关拉垮」。

---

## 阶段 2 · 内部验证独立性（在 TCGA 内，无需新数据）

### 2.1 leave-one-cancer-out (LOCO)

- 对 32 癌种的 M2/Myeloid 关联，逐癌种留出后重算荟萃。
- 观测：合并 ρ 是否稳定、CI 是否始终跨/不跨 0。
- 若删除某癌种后结果翻转，说明异质性由该癌种主导（已见 I²≈77%）。

### 2.2 预设验证癌种分离

- 在写稿时预先指定一组「训练癌种」与「验证癌种」。
- 只在训练集定阈值/模型，在验证癌种检验。
- 不强求必须显著，目的是证明「不是只在特定癌种上碰对」。

### 2.3 交叉验证 repeat（可选）

- 若内部有样本级 PSI，可做分层随机 split（seeded），在同一癌种内复现。
- 说明内部复现性，而非外部泛化。

---

## 阶段 3 · 真正的外部 isoform-level 验证（最有效、但需权限/成本）

### 3.1 获取转录本或 junction 级数据

候选（按可达性排序）：

| 数据源 | 级别 | 可达性 | 备注 |
|---|---|---|---|
| GDC 非 TCGA STAR junction（如 GLASS） | junction PSI | dbGaP 受控 | 需授权，最接近 isoform |
| CPTAC-3 splice-junction | junction PSI | 受控访问 | 同上 |
| 公开 long-read RNA-seq（ONT/PacBio） | transcript | 少见 | 若存在则最佳 |
| 本地合作队列 RNA-seq BAM | junction/transcript | 需合作 | 现实路径 |
| 已发表 isoform RT-PCR / 原位 | probe | 文献提取 | 作为独立证据而非 new data |

### 3.2 样本级 Junction-inclusion 验证

- 对 FL 关键外显子边界计算 PSI；
- 同一免疫评分；
- 样本级 Spearman；
- 校正 ZP3 总表达后重算（验证「比例增量超越总表达」）。

### 3.3 若获得 dbGaP/CPTAC 授权（记录为公司未来方向）

- 明文写入 Data availability 与 Limitations：
  「当前公开数据无法提供 isoform 级外部验证；已申请 dbGaP/CPTAC-3，作为后续独立验证」。
- 若不申请，则如实写「未获取」而非虚构。

---

## 阶段 4 · 论文定位与措辞（立即可以写）

### 4.1 把外部结果表述为「external gene-level non-replication」

推荐措辞（英文）：

> The gene-level ZP3–immune association was not reproduced in two small independent GBM cohorts. Because these datasets lacked transcript- or junction-level measurements, this analysis tested generalization of the gene-level biological premise rather than validation of the FL/RI proxy itself. The result therefore limits cross-cohort generalization but does not directly refute the isoform-specific association.

### 4.2 把 Limitations 改写成结构化的三层结论

1. **null 的来源**：n=32/24 低功效 + ZP3 低表达接近检测底 + GBM 异质性（I²≈77% 已有）；
2. **null 检验的对象**：gene-level 总表达，非 isoform proportion；
3. **对论文的意义**：强化「不能把 gene-level ZP3–免疫关联过度解读」的方法学结论，而非推翻内部 isoform 关联。

### 4.3 把 null 转化为「方法学正向论据」

- 强调：**总基因表达无法替代 isoform-level 测量**，这是一个有意义的正向发现；
- 可引用本组内部证据：FL/RI 比例与总表达解耦、组成控制后关联仍存在；
- 结论措辞定为：**「FL/RI 比例是一个具有内部多尺度支持、但外部基因级泛化受限的 bounded proxy」**。

---

## 阶段 5 · 冻结与审计（所有新分析必须走同一流程）

| # | 动作 | 产出 |
|---|---|---|
| 5.1 | 新增 `freeze_a3_external_null.py`（含诊断 + 功效 + 敏感性 + 留一） | `a3_external_null_diagnostics.csv` |
| 5.2 | 若做 LOCO，新增 `freeze_a3_loco.py` | `a3_loco_frozen.csv` |
| 5.3 | 更新 `freeze_a3_external.py` 自检，把「方向一致性」改为「如实记录 null」 | 冻结表字段注明 Analysis=diagnostic |
| 5.4 | 更新 manuscript v0.3.2：Results/Limitations/Data availability | 稿件修订 |
| 5.5 | 更新 `verify_manifest.py` 清单（新增文件） | 全量重跑 PASS |
| 5.6 | git 提交（conventional commits） | 新 commit hash |

---

## 优先级总结

| 优先级 | 做什么 | 时间/成本 | 对论文增益 |
|---|---|---|---|
| P0 | 阶段 1 诊断（检测质量+功效+敏感性） | 低，纯计算 | 高（直接解释 null） |
| P0 | 阶段 4 措辞改写 | 低 | 高（把 null 变方法学论据） |
| P1 | 阶段 2 LOCO 内部分离 | 中 | 中（证明内部稳健，非外部泛化） |
| P1 | 阶段 3.1 dbGaP/CPTAC 申请 | 中（流程） | 高（真正的 isoform 级验证） |
| P2 | 阶段 3.2 junction PSI 外部验证 | 高（需权限） | 最高（补齐最后缺口） |

> 核心纪律：**只做能如实冻结进仓库的分析**；null 状态下的每一条输出都写进
> Limitations；绝不为了「好看」篡改统计。