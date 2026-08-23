## 1. 总体评价

本文试图建立一个以 ZP3 全长转录本比例为核心的、具有明确边界的 bulk RNA-seq 关联性代理指标，并从肿瘤—正常差异、SpliceSeq 对照、通路分析、泛癌免疫关联、组成数据控制及外部重分析等多个层面进行评估。研究问题具有一定方法学价值，尤其是作者主动区分了生态层面一致性、样本层面一致性和平台依赖性。  
但目前稿件存在若干可能影响结论可信度的基础问题：核心样本和运行编号存在内部矛盾；肿瘤—正常比较及混合效应模型的批次、组织来源和混杂控制不足；“外部验证”“read-length effect”“internal transportability”等表述明显超过现有证据；此外，转录本比例的计算、组成数据处理、meta-analysis 和 kallisto 外部分析尚未达到可复现和可审计的程度。建议大修，且应优先重算关键分析，而不仅是文字降调。

---

## 2. Major issues（P0，不改可能被拒）

### Major 1：样本总数、分组和外部 FASTQ 运行数存在内部不一致

**问题—位置：**  
- Abstract、Results、Methods：总样本数为 19,131，其中 tumor 9,186、normal 7,792，合计仅 16,978，剩余 2,153 个样本被概括为“other”，但未说明其来源、类别及是否进入任何分析。  
- External isoform-level re-analysis：文中称“24 runs, SRR7050121–SRR7050184”，但该连续 accession 区间并非 24 个 run。  
- 文中同时使用“24 FASTQ files”“24 samples”“24 IDH-wild-type GBM”，但没有提供 sample–run 映射，也没有说明技术重复是否合并。

**具体修改建议：**  
1. 增加一张完整的样本流程表：TCGA、TARGET、GTEx 各自样本数，tumor/normal/other 定义，进入每个分析模块的样本数。  
2. 对 GSE113474 提供 `GSM–SRR–FASTQ` 映射表，明确：
   - 24 个生物学样本分别对应多少个 run；
   - 是否存在技术重复；
   - 技术重复如何合并；
   - 是否确实只分析 24 个样本，而不是一个 accession 区间中的全部 run。  
3. 在图、表和正文中统一 N；所有筛选和缺失值排除均应能由流程图重现。  
4. 若“SRR7050121–SRR7050184”只是一个较大 accession 范围中的部分 run，应列出实际下载的 accession，而不能用连续范围代替。

---

### Major 2：肿瘤—正常比较严重混杂于组织来源、癌种和平台，不能直接解释为“tumor-versus-normal switch”

**问题—位置：** Results“ZP3 transcript proportions shift between tumor and normal tissue”；Methods“Tumor-versus-normal comparison”。

TCGA/TARGET 肿瘤样本与 GTEx 正常样本在组织来源、取材条件、RNA 提取流程、测序平台、样本质量、个体构成和癌种分布上均可能存在系统差异。将 9,186 个 tumor 与 7,792 个 normal 直接进行 Mann–Whitney U 检验，得到极小 P 值并不能区分肿瘤生物学效应、组织组成差异和数据来源批次效应。对于 ZP3 这类组织特异性/异位表达基因，这一问题尤其严重。

此外，“tumor versus normal”是否在同一癌种内比较并不清楚；GTEx 正常组织与 TCGA adjacent normal 也不能简单视为同一类 normal。

**具体修改建议：**  
1. 将主分析改为癌种内、来源匹配的比较，至少报告：
   - TCGA tumor vs TCGA adjacent normal；
   - 对每个癌种分别分析；
   - TCGA/GTEx source-stratified 分析。  
2. 对肿瘤和正常样本加入数据来源、组织类型、癌种等协变量，或采用分层模型/批次校正后分析。  
3. 报告每个癌种的效应量、置信区间和方向一致性，而不只报告 pooled P 值。  
4. 将结论由“tumor-versus-normal switch”改为更谨慎的“在所分析数据集及其组织/平台构成下观察到的 tumor-associated proportion difference”，除非完成充分的匹配和敏感性分析。  
5. 明确“other”样本是否被排除，以及排除规则是否在分析前预设。

---

### Major 3：混合效应模型不足以支持“独立于总 ZP3 表达”的结论，且缺少关键模型诊断

**问题—位置：** Results“The FL proportion associates with immune features”；Methods“Mixed-effects models”。

当前模型似乎仅包括 FL proportion 或 RI proportion、癌种随机截距和免疫分数。至少存在以下问题：

1. 癌种随机截距不能自动消除测序来源、肿瘤纯度、组织学亚型、批次、免疫评分算法和临床特征等混杂。  
2. 免疫评分本身可能与癌种及样本来源高度相关，线性关系、方差齐性和随机效应结构均未验证。  
3. “调整总 ZP3 表达后仍显著”不能等价于“isoform usage carries information beyond total expression”。总表达与比例具有数学相关性，且 `total expression` 与比例共享同一组 TPM，可能产生严重共线性和解释歧义。  
4. 未报告 β 的量纲。比例从 0–1 变化一单位并不具有实际可解释性。  
5. 文中 Results 报告调整后 cytolytic P=未报告/0.10，而 Methods 中为 P=0.10；需统一。  
6. 仅报告 P 值，没有完整的 β、SE、95% CI、随机效应方差、边际/条件 R² 和模型诊断。

**具体修改建议：**  
1. 重新定义模型并明确公式，例如：
   \[
   \text{immune score}_{ij} = \beta_0+\beta_1\log(FL/RI)_{ij}+\beta_2\log(\text{total ZP3}+1)_{ij}+u_j+\epsilon_{ij}
   \]
   同时说明是否加入 cancer type、data source、purity、sex、batch、subtype 等协变量。  
2. 由于比例是有界变量，至少进行以下敏感性分析：
   - log(FL/RI)；
   - logit(FL proportion)；
   - beta regression 或稳健线性模型；
   - 将 FL/RI 绝对 TPM 与总 ZP3 分开建模。  
3. 报告 FL proportion 与 total ZP3 的相关性、VIF/条件数及模型残差诊断。  
4. 将 β 转换为有意义的效应量，例如 FL proportion 增加 0.1 或 FL/RI 比值翻倍时免疫分数的变化。  
5. 对癌种异质性至少加入 random slope，或报告固定效应模型与癌种分层结果；不要仅凭随机截距宣称跨癌种普遍性。  
6. “partial independence”“isoform usage, not merely total expression”应改为“在该模型和数据中，比例变量在加入总 ZP3 后仍携带统计学残差信息”，避免生物学独立性的暗示。

---

### Major 4：多重检验和主要终点定义不够清楚，当前“六个/七个显著”可能存在选择性叙述

**问题—位置：** Abstract、Results、Methods“Statistical inference and multiple testing”。

稿件包含大量分析和多个相关性体系：7 个转录本、7 个免疫特征、多个比例/比值、两类外部队列、GSEA 多个通路、SpliceSeq 多个事件和癌种分层 meta-analysis。虽然 Methods 提到 BH 校正，但没有清楚列出每个分析模块的检验总数、校正后的 q 值，也没有说明哪些终点是预先指定的、哪些是事后选择的。

尤其是：
- 主文中大量报告 nominal P；
- “六 of seven”被多次作为主要结果，但没有统一给出 q 值；
- M2 和 myeloid 被称为“headline immune features”，但这一选择是否预先设定不清楚；
- 外部分析包含多个特征和多种表达变换，不能只报告“all P>0.4”而不说明整个检验家庭。

**具体修改建议：**  
1. 增加“analysis plan / multiple-testing table”，列出每个模块：
   - 检验数量；
   - 校正家庭；
   - 校正方法；
   - 最终 q 值；
   - 主要/次要终点。  
2. 主文所有核心关联同时报告 β/ρ、95% CI、原始 P 和 q 值。  
3. 对“六个特征显著”改为基于预先定义阈值的正式结果，而不是按显著性数量进行概括。  
4. 对 compositional controls、LOCO、L2CO 和 held-out split 明确标记为敏感性/探索性分析，不应与 confirmatory primary analysis 混在一起。  
5. 对外部小样本分析重点报告效应量和置信区间，而非使用“null”作为二元结论。

---

### Major 5：癌种 meta-analysis 的统计方法与“transportability”表述不匹配

**问题—位置：** Results“cross-cancer meta-analysis”；Methods“Compositional controls and cross-cancer meta-analysis”。

文中使用每癌种 Spearman ρ，经 Fisher z 转换后以 N−3 加权，并采用 fixed-effect meta-analysis；但报告 I²=77%，提示明显异质性。固定效应 pooled estimate 在此情形下并不适合作为总体可迁移效应。更重要的是，癌种内相关性受到样本量、范围限制、癌种特异性免疫评分和技术差异影响，不能据此证明“consistent positive across cancers”或“internal transportability”。

LOCO/L2CO 分析并不能创造独立验证；held-out split 仍来自同一 TCGA 数据资源和同一处理流程。若验证集在事后选择，称为“prespecified”也需要提供预先注册或时间戳证据。

**具体修改建议：**  
1. 同时报告 random-effects meta-analysis、τ²、95% prediction interval 和各癌种效应量。  
2. 报告每个癌种的样本数、ρ、95% CI、方向和异质性来源。  
3. 不要把固定效应 pooled ρ 的精确 CI 解释为跨癌种普适性证据。  
4. 将“establish internal transportability”改为“在同一 TCGA 资源内的内部重采样/分癌种稳定性分析”。  
5. 说明 22/10 cancer split 的具体癌种、随机种子、确定时间以及是否在主要结果形成后才选择。  
6. 如果验证集 CI 由 pooled correlation 计算，应明确其方法；不应将同一数据源的 held-out 结果称为 external validation。

---

### Major 6：SpliceSeq 的“independent validation”及生态相关性被过度解释

**问题—位置：** Abstract、Results“The FL proportion tracks independent splice-event measurements”。

SpliceSeq 与 TCGA 样本存在数据来源或样本重叠，因此“independent”至少需要限定为“independent quantification resource”而不是 independent cohort。癌种层面 n=8 的均值相关性（ρ=0.95）容易受组织构成、癌种差异和均值聚合影响，属于 ecological concordance，不能证明样本级测量准确性。文中虽然提到 ecological correlation，但随后仍使用“reflects true splicing”“strongly correlated”进行较强推断。

此外，Methods 称有三个 AP events，但正文主要使用 AP1；应明确 AP1、AP2、AP3 与 FL/RI 转录本的对应关系，说明样本匹配、缺失值和每癌种样本数。

**具体修改建议：**  
1. 将“independent splice-event measurements”改为“independently generated/released SpliceSeq event measurements”，并明确样本是否重叠。  
2. 对 8 个癌种分别报告：
   - n；
   - FL mean/median；
   - AP PSI mean/median；
   - sample-level ρ；
   - 缺失率。  
3. 对生态相关性报告 leave-one-out 和 bootstrap，但同时强调 n=8 的不确定性；建议加入 permutation test 或 cluster-aware bootstrap。  
4. 结论改为“cohort-level concordance”，不能写成“true splicing validation”。  
5. 解释为什么 AP1 被作为主要事件，并预先定义 AP1 与 FL transcript 的生物学/结构对应关系。

---

### Major 7：GSEA 分组策略和排名方法不充分，存在组成、癌种和免疫浸润混杂

**问题—位置：** Results“FL-high and retained-intron-high samples differ in pathway enrichment”；Methods“GSEA”。

“FL-high versus retained-intron-high”没有定义阈值、是否互斥、是否按癌种内分组、每组样本数以及是否排除了低信号样本。若按照全体 GBM+LGG 排序，GSEA 很可能反映 GBM/LGG 的癌种差异、肿瘤纯度或细胞组成，而非 ZP3 isoform proportion。

此外，排名是“difference in expression (z-scores) between groups”，但没有说明是单基因 z-score、limma/t-statistic、Wilcoxon statistic 还是简单均值差。Hallmark 2020 版本、基因 ID 映射、重复基因处理和置换策略也未完整说明。

**具体修改建议：**  
1. 明确定义 FL-high/RI-high：
   - 阈值；
   - 组间是否互斥；
   - 每组 n；
   - 是否癌种内标准化/分层。  
2. 用 limma/DESeq2 等模型产生带有方差估计的全基因排名，而不是简单 z-score 差值。  
3. 加入癌种、纯度或免疫评分等协变量，或分别在 GBM 和 LGG 内进行 GSEA。  
4. 报告完整 pathway table，包括 NES、nominal P、FDR、leading-edge genes 和 gene-set size。  
5. 删除“FL-high enrichment is consistent with immune-context association”这类近似验证性措辞，改为“与相关性结果呈现方向一致的转录组模式”。

---

### Major 8：外部 kallisto 分析不能单独证明“read-length-dependent assignment effect”

**问题—位置：** Results“External isoform-level re-analysis”；Discussion；Methods“Kallisto”。

GSE113474 只有一个短读长平台和一个外部队列，同时还存在癌种/亚型、建库、测序平台、参考注释、定量软件（TCGA RSEM vs kallisto）等多个变化因素。因此 FL proportion 从 0.403 降至 0.012 不能被唯一归因于 read length。文中“demonstrates read-length/platform dependence”以及“most reads… therefore assignment favours RI”属于机制性解释，但没有模拟或基准实验支持。

特别需要注意：保留内含子转录本并不必然“包含所有外显子和内含子”且被 kallisto 自动偏向；转录本兼容性、有效长度校正、参考注释以及多重匹配读段均可能影响结果。

**具体修改建议：**  
1. 将结论降为：“该外部队列的 51-bp single-end/kallisto/Ensembl release 110 组合产生了与 TCGA RSEM 结果不同的比例，提示平台和定量流程依赖性；本研究无法区分读长、软件、注释和队列生物学因素。”  
2. 进行最基本的技术敏感性分析：
   - kallisto 与 Salmon/RSEM 的比较；
   - 不同 Ensembl release/自定义 transcriptome；
   - 改变 `-l/-s` 和检测阈值；
   - 是否加入 decoy sequence；
   - transcript-level posterior uncertainty 或 bootstrap。  
3. 进行模拟读段分析：从已知 FL/RI 混合比例的转录本生成 51-bp single-end 和 TCGA 类 paired-end reads，评估比例偏差。  
4. 检查区分 FL、RI 和 AP 的关键 junction 是否有支持读段，并报告每个样本的有效信息量，而不只报告 pseudo-alignment rate。  
5. 将“RI dominant”改为“在该量化配置下，RI transcript 获得较高的 estimated TPM”，避免将估计结果等同于真实生物学丰度。  
6. 不应把 24 个样本的外部重分析称为“external isoform-level validation”；更准确的名称应是“platform-specific technical re-analysis”。

---

### Major 9：转录本比例计算方法不够严谨，且正文与 Methods 可能不一致

**问题—位置：** Methods“Isoform proportion calculation”。

文中称：
- 先从 `log2(TPM+ε)` 矩阵反变换；
- 再以七个 ZP3 转录本 TPM 之和作为分母；
- 但前文又多次写“sum of TPM across all ZP3 transcripts”。

需要明确：
1. 是否原始输入就是 RSEM TPM，还是经过 log 变换后再反变换；
2. ε 的具体值；
3. 是否对零值进行了处理；
4. 分母究竟是七个指定转录本还是注释中的全部 ZP3 transcript；
5. 不同数据集是否使用完全相同的转录本版本和 transcript ID；
6. TPM 是否适合作为转录本比例的分子分母，特别是不同长度和 retained-intron transcript 的有效长度不同。

此外，`log2(TPM+ε)` 后反变换如果 ε 或四舍五入处理不当，可能改变低表达转录本的比例，尤其是外部 GSE113474 数据。

**具体修改建议：**  
1. 给出精确数学公式和伪代码。  
2. 明确七个 transcript 的完整列表、Ensembl release、canonical/retained-intron 注释来源和版本。  
3. 直接使用未变换 TPM 计算比例，并把“从 log 矩阵反变换”的处理作为敏感性分析。  
4. 对零值和低信号样本提供检测率、分母分布和比例不确定性。  
5. 分别报告：
   - seven-transcript denominator；
   - all annotated ZP3 transcripts denominator；
   - raw counts/estimated counts；
   - TPM-based estimates。  
6. 不要把该指标称为 PSI；建议全文统一为“ZP3 transcript fraction”或“FL transcript proportion”，仅在首次出现时说明其为 PSI-like，而非 splice-junction PSI。

---

### Major 10：外部“null”结果的解释过强，统计功效和结论边界需重写

**问题—位置：** Results“External cohorts”；Discussion limitations。

n=32 和 n=24 的外部队列确实不足以排除小到中等效应，但“all P>0.4”“null result”“most parsimonious explanation is limited power”仍然不能支持“null reflects limited power rather than measurement artifact”。零值比例为 0、IQR 合理和正对照显著，只能说明数据并非完全失效，不能证明 ZP3 关联不存在测量偏差或队列不匹配。

此外，两个外部队列均为 GBM，且可能具有不同 IDH 状态、治疗背景、样本来源、表达平台和免疫组成，不能作为对泛癌结果的直接否定或支持。

**具体修改建议：**  
1. 报告每个特征的效应量和 95% CI，而不只报告 P 值。  
2. 用“未观察到统计学显著关联，置信区间较宽，无法排除小/中等效应”替代“null”。  
3. 功效分析应说明效应假设、显著性水平、是否考虑多重检验以及计算方法。  
4. 将“not a measurement artifact”改为“未发现明显的零值或基本质量问题，但不能排除平台、批次、队列组成和样本量不足造成的差异”。  
5. 若可能，增加一个独立、样本量更大的 GBM 或其他癌种队列；否则应将外部结果作为限制，而不是“boundary on generalization”的确定性证据。

---

### Major 11：标题、摘要和结论仍然超过证据支持范围

**问题—位置：** Title、Abstract、Discussion、Conclusion。

当前标题“across bulk RNA-seq cohorts”容易让读者理解为该 proxy 已在多个独立队列中验证，但真正的主要证据来自 TCGA/TARGET/GTEx 内部资源；外部 GBM 基因级结果为未复制，外部 isoform-level 仅为技术重分析，且不可与 TCGA 直接比较。

“full-length transcript is proportionally elevated in tumors and in immune-featured samples”“isoform usage tracks immune context”“cross-cancer transportability supportive”等表述仍可能被理解为生物学泛化或机制证据。

**具体修改建议：**  
建议标题改为类似：

> **A platform-dependent ZP3 transcript-fraction analysis in TCGA/TARGET/GTEx reveals bounded associations with bulk immune-context scores**

或中文含义：

> **TCGA/TARGET/GTEx 中 ZP3 全长转录本比例与 bulk 免疫特征的有限关联及其平台依赖性评估**

摘要中应明确：
- 主要发现来自同一公共资源内的分析；
- 样本级 SpliceSeq 一致性有限；
- 外部基因级关联未复制；
- 外部短读长分析不能作为正式 isoform validation；
- 不应使用“validity”“transportability”“confirm”而不加限定。

---

### Major 12：可重复性声明目前不满足方法学期刊要求

**问题—位置：** Data availability、Code availability、Methods。

文中多次引用 `freeze_a3_*.py`、冻结结果、审计表和“project repository”，但目前没有公开 DOI、URL、commit hash、环境文件或完整命令。Code availability 使用“will be deposited prior to publication”，对方法学论文不足以支持审稿和复核。

还缺少：
- R/Python 版本；
- 包版本；
- kallisto index 的 hash；
- reference FASTA/GTF 的 checksum；
- 原始下载文件的 accession 和 checksum；
- 参数完整性；
- 样本排除记录；
- 免疫评分基因集及其版本；
- TCGA SpliceSeq 数据版本和下载日期。

**具体修改建议：**  
1. 在修回稿中提供可访问的公共仓库和 DOI，不要仅写“将发布”。  
2. 提供：
   - README；
   - 一键运行脚本或 workflow；
   - `environment.yml`/Docker/Singularity；
   - reference checksum；
   - kallisto index checksum；
   - 原始输入和冻结输出的 manifest；
   - 每一步的命令行参数。  
3. 将所有脚本、冻结表、图形源数据与正文编号一一对应。  
4. 对外部原始 FASTQ 提供实际下载日期，避免未来日期或无法核验的时间戳。  
5. 如果受版权或文件大小限制无法发布某些数据，应提供生成方式和公共 accession，而不是仅提供派生矩阵。

---

## 3. Minor issues（P1/P2）

### Minor 1：效应量和 P 值报告不完整

**问题—位置：** Fig. 1、Fig. 4、Supplementary Tables。

Fig. 1 仅给出 FL 的 effect r，未给 RI 和内部启动子转录本的效应量；Fig. 4 主要展示 P 值和 β，但缺少 SE、95% CI 和 q 值。  
**建议：** 所有主要结果统一报告 `N、effect estimate、95% CI、raw P、adjusted q`；图中可用点估计和 CI 替代大量指数形式 P 值。

---

### Minor 2：比例中位数与 TPM 结果之间的关系需解释

**问题—位置：** Fig. 5、Results。

“FL proportion median=0.012”“RI TPM median=4.24、FL TPM median=0.33”“median FL/RI ratio=0.031 in nine samples”不能直接相互换算，因为比例分母可能包括七个转录本，且 median ratio 不等于 ratio of medians。  
**建议：** 明确分母定义，增加 per-sample scatter/boxplot，并同时报告七个转录本总 TPM、FL/RI 检测率和比例的中位数/IQR。

---

### Minor 3：低信号过滤阈值需要客观依据

**问题—位置：** Compositional controls。

“dominant isoform proportion ≥0.5”被描述为低信号过滤，但该指标本身是相对量，不能直接代表 ZP3 的低表达；且过滤掉 1,609 个样本的依据不清楚。  
**建议：** 以总 ZP3 TPM、有效 reads、posterior uncertainty 或 transcript-level count 作为低信号标准，并将 ≥0.5 作为独立敏感性分析，而不是称作“low-signal filter”。

---

### Minor 4：log(FL/RI) 的数学表述不准确

**问题—位置：** Results、Methods。

文中称 `log(FL/RI)` 是“the simplest isometric log-ratio contrast”。对于两个组成部分，二元 log-ratio 与 ILR 仅差一个比例常数；但若七个转录本构成完整 composition，则 pairwise log-ratio 不等同于完整 ILR。  
**建议：** 改为“pairwise log-ratio contrast”；如使用 ILR，请给出完整正交基及公式。

---

### Minor 5：GSEA 的组名、方向和图例需统一

**问题—位置：** Fig. 3、Methods。

FL-high 与 RI-high 的 reciprocal pattern 可能只是同一轴的两端，不能视为两个独立结果。  
**建议：** 明确比较方向，例如“FL proportion-associated ranking”或“FL-high versus RI-high”；报告每个比较的 NES、FDR 和组别定义，避免重复叙述 reciprocal enrichment。

---

### Minor 6：免疫评分管线不够透明

**问题—位置：** Methods、Supplementary Table S3。

未列出七个免疫 feature 的具体基因集、算法、z-score consensus 的公式、是否在每个癌种内标准化、VSIR 缺失后的处理以及不同平台上的 gene-symbol 映射规则。  
**建议：** 将完整基因集、版本、缺失基因、标准化流程和代码作为补充材料；说明外部队列与 TCGA 是否使用完全相同的 gene set 和表达尺度。

---

### Minor 7：外部队列的引用和日期存在明显问题

**问题—位置：** External cohorts、References。

正文仍有 `[ref GSE77530]` 和 `[ref GSE113474]` 占位符。Methods/Data availability 的下载日期为 2026-08-18 至 2026-08-22，若与投稿日期不一致或属于未来日期，将严重影响可核验性。  
**建议：** 替换占位符，核对所有 accession、发布日期和实际下载日期；删除或更正未来日期。所有参考文献，尤其是 2026 年文献和与 ZP3 机制相关的引用，应逐条核验 DOI、PMID、卷页和是否已正式发表。

---

### Minor 8：参考文献与论证关联不够精准

**问题—位置：** Introduction、Discussion。

大量引用与 ZP3、isoform quantification 或 immune proxy 关系较弱，可能造成“引用堆砌”。部分文献似乎用于支持并未直接证明的机制性表述。  
**建议：** 精简背景引用，优先使用：
- ZP3 transcript biology；
- isoform quantification benchmarking；
- compositional data；
- bulk deconvolution limitations；
- read-length/platform effects；
- external validation methodology。  
每条引用应与具体论断直接对应。

---

### Minor 9：语言上应减少“confirm”“demonstrate”“establish”

**问题—位置：** 全文。

例如“compositional controls confirmed”“this analysis demonstrates read-length dependence”“establish internal transportability”“confirm the pooled estimate”。这些词在当前观察性、同源数据和小样本外部分析下过强。  
**建议：** 统一替换为“supported”“was consistent with”“suggested”“provided sensitivity evidence”，并在摘要、讨论和图例中保持一致。

---

### Minor 10：需要区分“detected”“quantified”“reliably quantified”

**问题—位置：** Fig. 5、Results、Discussion。

ZP3 TPM > 0 可能仅表示模型分配了非零估计量，并不等同于可靠检测。  
**建议：** 预先定义 detection threshold、minimum TPM/count、posterior uncertainty 或 bootstrap CV，并分别报告“non-zero estimate”“above threshold”“isoform confidently quantified”。

---

### Minor 11：外部单端 kallisto 参数需要说明其合理性

**问题—位置：** Methods。

`--single -l 51 -s 15` 的 fragment length 参数应与 GEO/SRA 的实际文库说明核对；如果 51 bp 是 read length 而非 fragment length，二者不能混用。  
**建议：** 明确 read length、fragment length、insert size 和 SD 的来源，并说明是否进行 adapter trimming、质量过滤和 poly-A 处理。

---

### Minor 12：标题和关键词中的“percent-spliced-in”可能误导

**问题—位置：** Title、Keywords、Abstract。

该指标不是 junction-defined PSI。  
**建议：** 在标题和关键词中优先使用“transcript proportion”“isoform fraction”或“PSI-like transcript fraction”，并在摘要第一处明确非标准 PSI。

---

### Minor 13：图表应展示癌种异质性，而不只是 pooled 结果

**问题—位置：** Fig. 4、Supplementary Table S2。

当前主图可能使读者看到总体正向效应，却无法判断多少癌种方向相反、哪些癌种样本量主导结果。  
**建议：** 增加按癌种分层的 forest plot 或 caterpillar plot，并展示 prediction interval、方向翻转比例及每癌种 n。

---

### Minor 14：伦理声明需核对 GEO 数据的具体使用条款

**问题—位置：** Ethics approval and consent。

公开数据不必然等同于所有衍生分析无需伦理审查，尤其是跨数据库匹配时。  
**建议：** 说明未访问受控临床数据，且根据所在机构政策获得豁免/无需审批；如无机构判定，应避免绝对化表述。

---

## 4. 最省力修订路径（按性价比排序的 7 条行动清单）

1. **先修复所有数据与编号矛盾。**  
   统一 19,131、9,186、7,792、7,577、694、32、24 等所有 N；提供样本流程表和 GSE113474 的 GSM–SRR 映射。

2. **把核心结论全面降调并修正标题/摘要。**  
   将“validated”“demonstrates”“establishes transportability”“null”分别改为“proxy concordance”“technical discrepancy”“internal stability analysis”“未观察到显著关联”。明确外部 kallisto 结果不是正式 isoform validation。

3. **补充完整转录本比例定义并重新计算一套主分析。**  
   明确七个转录本、参考版本、ε、分母和零值处理；优先使用原始 TPM 计算比例，并将 log-ratio 作为主要敏感性分析，而不是事后补充。

4. **重做或至少完整诊断混合效应模型。**  
   加入总 ZP3、癌种及可获得的来源/纯度协变量；报告 β、SE、95% CI、q 值、VIF、随机效应和模型诊断；将“独立于总表达”改为“条件关联”。

5. **重新处理 tumor–normal 分析。**  
   至少增加癌种内、TCGA-only 和 source-stratified 结果；如果无法完成充分匹配，将 tumor–normal 结果明确标为跨资源描述性差异，而不是生物学 switch。

6. **把 meta-analysis 和 GSEA 改成可审计的补充分析。**  
   meta-analysis 增加 random-effects、prediction interval 和每癌种结果；GSEA 明确分组规则、排名统计量、每组 n 和癌种控制。

7. **立即公开完整复现材料。**  
   提供 DOI、commit hash、环境文件、reference/index checksum、命令、输入输出 manifest、冻结表和图源数据；同时核对所有参考文献、下载日期和 `[ref GSE...]` 占位符。

完成上述 7 项后，稿件可以较为可信地定位为：**一个展示 ZP3 transcript fraction 在特定 TCGA/TARGET/GTEx 资源中与 bulk immune scores 存在有限关联、且该关联受数据组成和测序平台影响的分析框架**，而不是已经验证的跨队列 isoform biomarker。