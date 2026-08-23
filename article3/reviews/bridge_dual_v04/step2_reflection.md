## 1. 总体判断

这份审稿意见总体上是**方向正确、科学立场基本公允、对文章最核心的风险判断准确**的一份意见。它抓住了真正的问题：这不是一篇 mechanistic validation 或 clinical biomarker-validation 研究，而是一个以 bulk RNA-seq 为基础、跨队列/跨癌种、具有 platform dependence 的 association-level transcript-usage study。审稿意见要求作者把“marker”“proxy”“immune-context association”“TREM2-associated biology”等概念分层，是合理且必要的。

不过，这份意见也有几个明显问题：

1. **把不同严重程度的问题都标成了 Major。**  
   例如 title shortening、reference pruning、workflow schematic、precision formatting，不应与 FDR family、外部复现定位、数值一致性放在同一层级。

2. **部分建议默认数据和分析管线仍然可获得。**  
   例如 alternative annotation、不同 kallisto 参数、random slopes、Cox Schoenfeld 检验等，未必都能在当前数据条件下完成。它们更适合写成“若要保留相应强结论，则建议补充”，而不是绝对必做。

3. **对若干未发生的分析存在条件性批评。**  
   关于 Cox proportional-hazards 和 grade encoding，如果稿件根本没有 survival analysis，这一段主要是“防止未来过度扩展”，而不是当前 blocker。将其作为现有统计缺陷会有些过度。

4. **对“marker”语言的警惕是对的，但有削弱文章贡献的风险。**  
   文章不应自称 validated biomarker，但也不必把结果降格为纯描述。若数据确实显示 tumor–normal isoform shift、cohort-level concordance 和一定 sample-level association，可以保留“bounded proxy”或“association-level descriptor”这一积极但克制的定位。

总体而言，我会评价为：**一份高质量但偏“防御性过度”的审稿意见；核心诊断准确，优先级需要重排，若逐条机械执行可能导致不必要的降格和分析扩张。**

---

# 2. 逐条反思

## A. 对审稿意见八大块的反思

### 1）Study design & logical framing

- **是否确属 blocker：部分是。**  
  如果摘要、标题、结论把 ZP3 FL proportion 写成 validated biomarker、myeloid-cell marker 或 mechanistic indicator，那么这是实质性问题；如果只是个别段落措辞偏强，则不是实验层面的 blocker，而是 narrative correction。

- **是否应该采纳：应该采纳。**

- **是否过度严格/可商榷：总体不过度，但“bystander marker”不必完全删除。**  
  只要明确定义为 non-causal, context-associated proxy，并说明不能推断 cell-of-origin、protein expression 或 function，就可以保留这一概念。完全删除可能损失文章的概念贡献。

- **minimal viable fix：**
  1. 统一全文主定位：  
     **“association-level, platform-dependent transcript-usage proxy associated with myeloid-enriched immune states in bulk glioma datasets.”**
  2. 首次出现时定义 marker/proxy 的含义和边界。
  3. 删除或改写“validated biomarker”“establishes”“mediates”等强词。
  4. title 不必立刻大幅重写，但应避免“isoform-resolved proxy for immune context”给人已完成 biomarker validation 的印象。

**判断：P0/P1 之间；对结论性文章而言接近 P0。**

---

### 2）Statistical methodology

#### A. Cell-level vs patient-level inference

- **是否确属 blocker：是概念性 blocker，但不是数据分析 blocker。**  
  如果把 bulk association 写成“ZP3 与 TREM2+ myeloid cells 共表达”或“ZP3 来自髓系细胞”，则会造成错误推断。

- **是否应该采纳：必须采纳。**

- **是否过度严格：不过度。**  
  这是 bulk transcriptomics 最基本的 ecological inference 限制。

- **minimal viable fix：**
  - 在 Methods、Results 首次描述 immune features 时写明：  
    “all associations are sample-level and derived from bulk RNA-seq.”
  - 在 Discussion 明确：不能识别 cell-of-origin，也不能证明同一细胞内 co-expression。
  - 不必因为没有 single-cell/spatial data 就否定 immune-context association；只需把它限定为 bulk-sample context association。

**判断：P0。**

#### B. FDR handling

- **是否确属 blocker：若正文确实无法追溯测试 family，是统计报告 blocker。**

- **是否应该采纳：必须采纳。**

- **是否过度严格：要求 family 定义合理，但要求所有结果都重新大规模校正可能过严。**  
  不同分析模块可以有不同 predefined family；没有必要把 tumor-normal、GSEA、pan-cancer、external cohort 全部强行并入一个总 family。

- **minimal viable fix：**
  - 给每个分析模块定义 family，例如：
    - tumor-normal transcript comparisons；
    - immune-feature associations；
    - pan-cancer analyses；
    - GSEA pathways；
    - external validation tests。
  - Supplementary tables 增加 raw P 和 FDR q。
  - 正文多重比较结果优先报告 q-value。
  - 若无法重算，至少明确哪些 P 是 nominal、哪些已 FDR-adjusted，并避免将 raw P 写成确证性结果。

**判断：P0/P1；比原审稿意见中大多数项目更重要。**

#### C. Cox proportional-hazards assumptions & ordinal grade encoding

- **是否确属 blocker：当前稿件没有 survival analysis 时，不是 blocker。**

- **是否应该采纳：有条件采纳。**  
  如果 manuscript 没有 survival claim，不需要新增 Cox analysis，也不应因未做 Schoenfeld test 而被判统计不合格。

- **是否过度严格/可商榷：原意见在这里有些“预防性过度审查”。**  
  它正确地指出 grade 不应未经论证直接作为线性/ordinal covariate，但这只在实际使用 Cox 或 grade modeling 时成立。

- **minimal viable fix：**
  - 若没有 survival analysis：删除 survival-adjacent wording，并在局限性中写明“survival prediction was not evaluated”。
  - 若已有 survival analysis：再补 PH assumption、grade categorical encoding 和关键分层/协变量说明。
  - 不建议为满足审稿意见而新增一套 survival analysis。

**判断：P2；除非文中确有 survival claim。**

#### D. Pan-cancer descriptive vs inferential claims

- **是否确属 blocker：不是分析 blocker，但属于重要措辞问题。**

- **是否应该采纳：应该采纳。**

- **是否过度严格：对“pan-cancer marker”的警惕合理；但不必把 pan-cancer analysis 降为毫无推断意义。**  
  可以保留跨癌种的 breadth，同时正面报告 heterogeneity，例如 I²≈77%。

- **minimal viable fix：**
  - 把“pan-cancer validation/generalization”改为：
    - “cross-cancer descriptive pattern”
    - “heterogeneous cross-cancer association”
    - “broad but non-uniform association”
  - 图中同时展示 effect size、方向和 heterogeneity。
  - 不需要重新做复杂 meta-analysis，除非原文将 pooled effect 当作强 generalization 证据。

**判断：P1。**

#### E. Threshold sensitivity in single-cell/short-read re-analysis

- **是否确属 blocker：取决于外部 re-analysis 在稿件中承担什么作用。**  
  如果它被称为 external validation，则是重要问题；如果只是说明 platform/read-length dependence，则不是 blocker。

- **是否应该采纳：方向应采纳，分析规模可酌情。**

- **是否过度严格：要求多个 annotation、参数和 benchmark 有些过重，尤其是短 reads 与原始数据/环境可能无法完全重建时。**

- **minimal viable fix：**
  - 明确“detected”与“reliably quantified”的区别。
  - 报告 FL 和 RI 同时超过最低 abundance threshold 的样本数。
  - 至少做一个 detection threshold sensitivity 或在补充材料中提供 threshold-dependent summary。
  - 若无法重算，使用更谨慎的表述：  
    **“consistent with read-length- and reference-dependent assignment effects”**，而不是 “demonstrates systematic assignment bias”。
  - 补充 annotation version、effective-length correction 和过滤阈值信息。

**判断：P1；若只是探索性分析，则 P2。**

---

### 3）Results reporting

#### Problem 1：Numerical consistency

- **是否确属 blocker：是编辑和可信度层面的 blocker。**

- **是否应该采纳：必须采纳。**

- **是否过度严格：不过度。**  
  Abstract、Results、figure legends、supplementary tables 中的 N、effect direction、P/q、测试对象必须一致。

- **minimal viable fix：**
  - 建立一个“canonical statistics table”。
  - 逐一核对：
    - sample size；
    - tested contrast；
    - effect direction；
    - raw P vs FDR q；
    - event 数量；
    - external cohort 的 exact N 和 power statement。
  - 不一定需要重新做分析，先做结果审计即可。

**判断：P0。**

#### Problem 2：Precision and overprecision

- **是否确属 blocker：不是。**

- **是否应该采纳：建议采纳。**

- **是否过度严格：略显吹毛求疵。**  
  小数位本身不会改变科学结论，除非结果极其接近 threshold。

- **minimal viable fix：**
  - correlations 报 2 位；
  - proportions 依据量纲统一 2–3 位；
  - P-value 用统一 scientific notation 或 2–3 significant digits；
  - 避免展示不具有稳定性的过度精确数值。

**判断：P2。**

#### Problem 3：Null/attenuated results

- **是否确属 blocker：对“external validation”叙事是 P1；对整篇文章不是 blocker。**

- **是否应该采纳：应该采纳。**

- **是否过度严格：不过度，但应避免把 null result 写成决定性否定。**  
  External gene-level null 可以削弱 transferability claim，但不必否定 isoform-level internal evidence。

- **minimal viable fix：**
  - 明确写出：external gene-level replication was null。
  - 说明这使“robust cross-cohort biomarker”结论不成立。
  - 将 isoform-level re-analysis定位为 exploratory/platform-specific evidence，而非真正 external validation。
  - Power limitation 可以保留，但不能成为“救回 null result”的唯一解释。

**判断：P1。**

---

### 4）Discussion restraint & attribution

- **是否确属 blocker：如果正文把文献机制写成本文数据已证实，则是重要 blocker；否则是 P1。**

- **是否应该采纳：必须采纳。**

- **是否过度严格：总体合理，但不必把所有机制相关术语全部删除。**  
  “TREM2-associated”可以作为 hypothesis-generating context；不能写成“ZP3 drives TREM2 biology”。

- **minimal viable fix：**
  将内容分成三层：
  1. **Observed in this study**：tumor-normal shift、sample-level association、heterogeneous immune association；
  2. **Literature-supported context**：GPX4–ZP3、myeloid biology、TREM2 background；
  3. **Speculation/future hypothesis**：cell-type specificity、functional consequence、membrane accessibility。
  
  推荐使用：
  - “consistent with”
  - “may reflect”
  - “hypothesis-generating”
  - “requires direct validation”

- **GPX4–ZP3 的处理：**  
  原审稿意见要求明确区分 mechanistic study 与 transcript-usage/proxy study，这一点应采纳。

**判断：P1，若有明确因果措辞则升至 P0。**

---

### 5）Figure–text–table consistency

- **是否确属 blocker：多数不是 blocker，但属于重要的可读性和可审计性问题。**

- **是否应该采纳：应采纳核心部分。**

- **是否过度严格：要求新增 workflow schematic 是可选项，不应作为硬性要求。**

- **minimal viable fix：**
  - Figure 4 legend 明确到底哪些 immune features 进行了 meta-analysis，哪些只是单队列结果。
  - Supplementary Table S2 和正文使用相同分析名称。
  - Figure 5 明确它展示的是 read-length/platform compatibility evidence，而非完整 benchmark。
  - 将 external null、diagnostics、LOCO/L2CO 等分开或在表头分区。
  - 不必新增图，只需修正 legend 和 Results 叙述即可。

**判断：P1/P2；legend 与表格错配属于 P1，workflow schematic 属于 P2。**

---

### 6）Language, structure, and formatting

- **是否确属 blocker：不是科学 blocker。**

- **是否应该采纳：应适度采纳。**

- **是否过度严格：把“过长”“reference list 太长”标为 Major，明显偏重。**  
  只要文章逻辑清楚，长 reference list 本身不会构成拒稿理由。真正的问题是重复性和结构负担。

- **minimal viable fix：**
  - 压缩重复出现的 proxy boundary。
  - 将方法细节移入 Supplementary Methods。
  - Discussion 围绕 4–5 个核心结果重组。
  - title 和 abstract 做适度精简。
  - reference audit 可以做，但不需要大规模删文献以“证明”概念清晰。

**判断：P2。**

---

### 7）Data/code availability & ethics statements

- **是否确属 blocker：对于最终发表可能是合规 blocker；对于科学结论不是。**

- **是否应该采纳：应采纳，但区分必需项和理想项。**

- **是否过度严格：要求现在就有 DOI/repository 很理想，但部分期刊允许 acceptance 后上传。**

- **minimal viable fix：**
  - 给出 repository URL、accession、代码版本和关键环境信息。
  - 明确哪些数据是 public、哪些是 controlled-access。
  - 提供 frozen result tables 或说明如何从 public inputs 重建。
  - ethics 改为：  
    **“This study used only public, de-identified data and did not require new ethics approval.”**
  - AI disclosure 按目标期刊 policy 调整；无需将 AI-use 与 scientific validity 混为一谈。

**判断：P1/P2；目标期刊有硬性数据政策时为 P0。**

---

### 8）Narrative positioning versus the cited Cell 2026 GPX4–ZP3 study

- **是否确属 blocker：不是独立 blocker，但对 novelty 和 attribution 很重要。**

- **是否应该采纳：应该采纳。**

- **是否过度严格：不过度。**

- **minimal viable fix：**
  在 Introduction/Discussion 各加一处清晰区分：
  - Cell study：mechanistic/pathway-oriented；
  - 本文：transcript-usage、proxy evaluation、cross-cohort association、platform dependence。
  
  明确本文不提供：
  - GPX4–ZP3 mechanism；
  - protein isoform validation；
  - spatial localization；
  - immune modulation function。

**判断：P1。**

---

## B. 对“7条 Major revision recommendations”的反思

原审稿意见中的 Priority list 实际上列出的是 **5 条 Highest priority + 4 条 Medium + 3 条 Lower priority**，并不是严格意义上的“7条 Major”。因此这里按最接近“主要修改事项”的前 7 项重排如下：

| 原 Priority 项目 | Blocker？ | 是否采纳 | 评价与 minimal viable fix |
|---|---:|---|---|
| 1. Downgrade causal/biomarker language | 对主结论而言是 | 必须 | 全文统一为 bounded association-level proxy；不要完全删除贡献 |
| 2. Clarify statistical families/FDR | 是统计报告 blocker | 必须 | 按分析模块定义 family，补 raw P/q/effect size/N |
| 3. Tighten external validation claims | 对验证叙事是 | 必须 | 将 external gene-level null 与 internal isoform analysis 分开；避免“validated” |
| 4. Resolve figure/text/table alignment | 是可信度问题 | 必须 | 做一致性审计、改 legends 和 supplementary table |
| 5. Add bulk-data limitations | 概念性 blocker | 必须 | 明确 bulk sample-level、cell-of-origin unknown、no co-expression proof |
| 6. Shorten title and Abstract | 否 | 建议 | 精简 title/abstract，但不必彻底压缩到牺牲信息量 |
| 7. Improve Methods transparency | 否，除非分析无法复现 | 建议 | 补 mixed model/GSEA/kallisto 关键参数；无法重算时以 limitation 替代 |

原列表中的后续项目也应保留，但优先级更低：

- null external gene-level replication：P1；
- reference trimming：P2；
- language/terminology：P2；
- workflow schematic：P2；
- precision formatting：P2。

---

## C. 对“7条 Minor revision recommendations”的反思

原文实际上列出 A–H 共 **8 条 Minor**。其中 A–G 可作为七条主要 minor 项目，H 另行处理。

| Minor 项目 | Blocker？ | 是否采纳 | 更务实的做法 |
|---|---:|---|---|
| A. Abstract wording | 否 | 建议 | “validated”改为 “evaluated/tested/benchmarked”；说明 platform-specific |
| B. Terminology consistency | 否 | 应采纳 | FL/RI/proxy/ecological association 统一，减少概念漂移 |
| C. SpliceSeq terminology | 否，但若事件对应关系不清会影响可解释性 | 应采纳 | 增加一小段 transcript/event mapping，最好配简图或补图 |
| D. Mixed-effects details | 否 | 视数据和模型而定 | 补 random intercept、standardization、residual checks；random slopes 不必强制 |
| E. GSEA details | 否 | 应采纳 | 报 ranking metric、gene-set size、permutation、tie handling |
| F. External re-analysis wording | 否 | 应采纳 | “consistent with/suggests/compatible with”，避免“most plausible explanation” |
| G. Supplementary tables | 否 | 建议 | 补 raw P、q、effect size、CI、exact N；CI 不可得时注明 |
| H. Reference list | 否 | 可选择性采纳 | 删除明显无关或重复引用，不必为压缩而机械删减 |

这里 D 项需要特别注意：  
“random slopes were considered”“residuals were checked”是透明度建议，不等同于必须重新拟合所有模型。若模型已经完成且结论不依赖复杂 random-effects structure，可以在 Methods 中说明当前 specification 的理由。

---

# 3. 优先级重排

## P0：不改可能导致拒稿或结论不被接受

### P0-1. 主结论和标题的 claim calibration

必须把“validated marker”“mechanistic marker”“cell-type marker”降为：

- association-level proxy；
- bulk-sample descriptor；
- platform-dependent transcript-usage correlate；
- myeloid-enriched immune-state-associated signal。

**理由：**这是文章 validity 的核心，而不是单纯语言问题。

### P0-2. FDR family 和统计报告可追溯性

必须让读者知道：

- 每一组测试的 family 是什么；
- P 是 nominal 还是 adjusted；
- q-value 如何得到；
- N 是否因过滤而变化。

**理由：**如果多重比较框架不清，所有“显著关联”都可能被质疑。

### P0-3. 数值、N、图表和 correction status 一致

Abstract、Results、legends、Supplementary tables 的关键数字必须统一。

**理由：**这直接影响编辑和审稿人对分析可靠性的判断。

### P0-4. 明确 bulk RNA-seq 的推断边界

不能让 sample-level association 被读成 cell-level co-expression 或 myeloid cell-of-origin 证据。

**理由：**这是 ZP3–TREM2/myeloid interpretation 最容易越界的地方。

### P0-5. 外部验证的定位

外部 gene-level null result 必须正面呈现；external isoform-level re-analysis 若存在 platform dependence，就不能称为一般意义上的 independent validation。

**理由：**“validation”是当前稿件最容易被攻击的词。

---

## P1：强烈建议修改，但通常不要求新实验

### P1-1. Pan-cancer 语言降级为 heterogeneous cross-cancer association

不必删除 32 cancer types 结果，但要避免“generalizes across cancers”的无条件表述。

### P1-2. TREM2、GPX4–ZP3 和其他机制背景分层

明确哪些是本文观察、哪些是文献背景、哪些是未来假说。

### P1-3. Short-read/external isoform analysis 增加最低限度 sensitivity information

优先补充 threshold、annotation version、共同检测样本数；没有条件时至少加 caveat。

### P1-4. Figure legend 和 supplementary organization 修正

尤其是 Figure 4 的七个 immune features 与 S2 只分析部分 features 的潜在错配。

### P1-5. Data/code availability 达到可复现标准

给出 accession、repository、版本和 frozen tables。若数据受控访问，说明获取路径。

### P1-6. Null result 的 Discussion 更对称

可以讨论 power，但必须先承认 null result 对 transferability claim 的削弱。

---

## P2：可改可不改，主要是编辑层面

- title 和 abstract 的进一步压缩；
- precision formatting；
- workflow schematic；
- reference list 大规模裁剪；
- random slopes、alternative annotation 等额外 sensitivity analyses；
- survival analysis/Cox Schoenfeld，如果原稿没有 survival claim；
- 逐项删除所有“marker”词，而不是重新定义其边界。

---

# 4. 对作者最省力的修订路径

建议不要把修订做成“重新分析论文”，而应采取以下最小可行路线。

## Step 1：先重写一句主结论

建议全稿围绕下面这句话组织：

> **ZP3 full-length transcript proportion is a bulk-sample, association-level and platform-dependent transcript-usage proxy that shows heterogeneous concordance with myeloid-enriched immune states in glioma datasets.**

这句话既保留贡献，也主动限定 inference boundary。

---

## Step 2：进行一次全稿 terminology pass

将下列表达统一：

- “validated marker” → “evaluated proxy” / “association-level proxy”
- “immune marker” → “immune-context-associated signal”
- “myeloid marker” → “myeloid-enriched context-associated feature”
- “TREM2-associated biology” → “association with TREM2-related/myeloid-enriched signatures”
- “systematic assignment bias” → “platform- and read-length-dependent assignment effects”
- “establishes” → “supports” / “is consistent with” / “suggests”

但不建议把所有“marker”全部删除。若使用，定义为 **non-causal context-associated proxy**。

---

## Step 3：补一个简短的统计说明，而不是全面重做

Methods 中新增一个“Statistical inference and multiple testing”小节，说明：

- 单位是 patient/sample；
- 所有 immune scores 来自 bulk RNA-seq derived signatures/deconvolution；
- 每类分析的 multiple-testing family；
- P 与 q 的区别；
- exact N 根据可用数据变化。

如果原始数据仍可获得，补充 q-value 和 effect size；如果不可得，至少不要把 nominal P 当作 confirmatory evidence。

---

## Step 4：用图注和局限性声明解决 cell-of-origin 问题

可以直接在 Figure 4/5 legend 或 Discussion 加：

> “Associations are computed at the bulk-sample level and do not establish cellular co-expression, cell-of-origin, spatial localization, or causality.”

这比新增 single-cell 或 spatial experiment 更符合当前研究设计。

---

## Step 5：把 external result 重新命名和分区

将 external analysis 分成：

1. **External gene-level assessment: null or attenuated replication**
2. **Platform-specific isoform re-analysis: exploratory compatibility assessment**

不要把后者统称为 external validation。图题和表题也要一致。

---

## Step 6：针对外部 kallisto 分析，优先做“轻量敏感性分析”

如果数据和环境还在，最值得补的不是全面 benchmark，而是：

- 两到三个 detection threshold；
- FL 和 RI 同时被检测的样本数；
- annotation/version；
- mapping/quantification quality summary。

如果无法重算，可以用以下方式替代：

> “Because the external short-read re-analysis was performed under a fixed historical annotation and limited read-length setting, its isoform proportions should be interpreted as platform-specific estimates rather than definitive measurements.”

这可以作为文字+局限性声明，不必硬补多个 pipeline。

---

## Step 7：做一次“canonical numbers audit”

建立一个内部表格，列出：

- N；
- comparison；
- effect size；
- rho；
- raw P；
- q；
- pathway/NES；
- exact cohort；
- figure/table location。

这通常比重新跑分析更能消除审稿人疑虑。

---

## 关于审稿意见提到的几个“重算”项目

### 1. Pseudobulk

如果原研究本身不是 single-cell cell-level analysis，而是 bulk RNA-seq association study，通常**不需要为了迎合审稿意见新增 pseudobulk**。可用 bulk sample-level limitation 解决。

### 2. Mixed-effects model

如果原模型和研究问题匹配，且数据不可得：

- 在 Methods 中说明 random intercept 的设定；
- 说明 sample size imbalance 和 standardization 的处理；
- 在 Discussion 中承认未系统比较 random slopes 或 alternative covariance structures。

不必把“considered random slopes”变成硬性重算要求。

### 3. Cox Schoenfeld

如果没有 survival analysis，最省力且最正确的方案是：

- 删除 survival-adjacent language；
- 说明 survival prediction/prognostic validation 不在本研究范围内。

不要为了满足意见而加入一项新的、与主线不一致的 survival analysis。

### 4. Alternative kallisto/annotation benchmarking

若无原始 FASTQ、旧版 reference 或运行环境：

- 不要声称完成了 robustness benchmark；
- 将结果定义为 exploratory platform-specific re-analysis；
- 增加 threshold、read length、annotation compatibility 的 limitation；
- 图注中加 caveat。

### 5. Confidence intervals

若 effect size 的 CI 能通过现有 summary data 计算，建议补；若不能，不要伪造精确 CI。可报告 effect size、P/q 和 exact N，并明确 CI 未计算的原因。

---

# 5. 潜在风险提醒

如果作者完全照单全收原审稿意见，存在明显的**过度削弱文章**风险。

## 风险一：把“bounded proxy”写成“没有生物学意义”

文章可以同时做到：

- 不声称 ZP3 是 validated biomarker；
- 不声称 ZP3 来自髓系细胞；
- 不声称 TREM2 是因果机制；
- 但仍然主张：ZP3 FL proportion 是一个有重复性但受平台约束的 transcript-usage descriptor，与特定 immune-context signatures 存在跨队列、异质但可解释的关联。

这不是“自我否定”，而是合理的 scientific positioning。

## 风险二：过度强调 external null，掩盖内部结果

External gene-level null result 应当明确报告，但不能写成“因此 ZP3 与 immune context 没有关系”。更准确的表达是：

- gene-level transferability was not demonstrated；
- isoform-level association remains suggestive；
- platform/read-length dependence may limit replication；
- this supports bounded rather than universal biomarker claims.

## 风险三：为补“严谨性”而加入不必要分析

Cox、pseudobulk、complex mixed models、多个 quantification pipelines，如果不是研究设计原本的一部分，可能会：

- 引入新的 multiple-testing burden；
- 造成方法堆积；
- 偏离 ZP3 transcript-usage 主线；
- 给审稿人更多可攻击的技术细节。

## 风险四：把“myeloid-enriched”误写成“myeloid-derived”

建议坚持以下区分：

- **myeloid-enriched context**：可以使用，前提是 bulk signature/deconvolution 支持；
- **myeloid-derived ZP3**：当前证据不足；
- **TREM2+ macrophage expression of ZP3**：需要 single-cell/spatial 或 protein-level co-localization；
- **ZP3 drives immunosuppression**：需要 functional evidence。

## 把握平衡的原则

最合适的定位不是“validated immune marker”，也不是“仅仅是统计相关”，而是：

> **a technically resolved, cohort-level transcript-usage signal that provides a testable proxy for immune-context association, with explicit limits in cell-of-origin, mechanism, and platform transferability.**

这个表述保留了：

- transcript-usage resolution；
- multi-cohort comparison；
- immune-context relevance；
- platform dependence；
- 可检验的后续假说。

---

# 6. 一句话总结

若我是编辑，我会将这份审稿意见定性为**核心判断准确、但把若干可选分析和编辑性问题过度升级为 Major revision 的偏严格意见**；编辑决定建议为：**Major revision，但以 claim calibration、FDR/数值可追溯性、外部验证定位和 bulk-data inference boundary 为必要条件，不要求新增 survival、single-cell/spatial 或全面重做分析。**