Below is a candid reviewer-style critique, prioritized and section-specific. Overall, the manuscript has a strong data-processing/comparative bioinformatics effort and an unusually careful instinct to state boundaries. However, the current narrative still overreaches relative to the evidence in several places, and the statistical reporting is not yet fully aligned with the inferential claims. In particular, the paper should be reframed as an association-level, platform-dependent transcript-usage study with limited sample-level validation, not as a “marker” study in the biomarker-validation sense.

---

# Overall assessment

**Major concern:** the manuscript currently oscillates between three different claims:

1. **ZP3 FL proportion is a biologically meaningful isoform-usage descriptor in TCGA/GTEx.**
2. **ZP3 FL proportion is a proxy for alternative-promoter/splice-event behavior at the cohort level.**
3. **ZP3 FL proportion is a myeloid/immune-context marker, potentially of TREM2-associated biology.**

The evidence is strongest for (1), moderate for a constrained version of (2), and only weak-to-hypothesis-generating for (3). The manuscript will improve substantially if it **downgrades the “marker” language** and clearly distinguishes:

- **descriptive association**
- **ecological concordance**
- **sample-level concordance**
- **platform-transferability**
- **biological mechanism**

Right now, some sections imply a level of validation that is not justified by the data.

---

# Major revision recommendations

## 1) Study design & logical framing
### Major
**Section(s): Abstract, Introduction, Discussion, Conclusion**

### Problem
The “context/bystander marker” conclusion is currently **too strong and too teleological** for the evidence presented. The study is association-level and multi-cohort, but it lacks:
- protein validation
- spatial validation
- single-cell co-localization
- perturbation/functional evidence
- direct evidence that FL ZP3 arises from myeloid cells rather than tumor cells, stromal cells, or bulk composition shifts

The manuscript repeatedly uses “marker,” “immune-context association,” “proxy marker,” and “boundary” correctly in some places, but elsewhere still approaches implied biomarker validation. The title also implies a more refined mechanistic status than the data support.

### Actionable fix
Reframe the central conclusion as:

- **“ZP3 FL transcript proportion is an association-level, platform-dependent isoform-usage proxy linked to myeloid-enriched immune states in bulk glioma datasets.”**

Avoid “bystander marker” unless you explicitly define it as a **non-causal, context-associated proxy**. If you want to keep that concept, you need one paragraph in the Discussion explaining why the term is used and what it does **not** imply.

### Also needed
The title is too long and conceptually overloaded. Consider simplifying to something like:

- **“ZP3 full-length transcript proportion is an association-level proxy for immune context in glioma across bulk RNA-seq cohorts”**

That would be more honest and more readable.

---

## 2) Statistical methodology
### Major
**Section(s): Materials and methods; Results; Supplementary Tables S1–S4**

This is the most important area for revision.

### A. Cell-level vs patient-level inference
#### Problem
You infer biological associations from bulk transcript proportions using sample-level statistics, but the manuscript sometimes slips into language suggestive of cell-state inference. The immune feature scores are deconvolution-derived or signature-derived bulk features, not cell counts. Likewise, isoform proportion is bulk transcript usage, not necessarily tumor-cell intrinsic usage.

### Actionable fix
State explicitly, at first use and in every relevant result/discussion paragraph:
- these are **bulk-sample associations**
- they do **not** identify the expressing cell type
- they do **not** prove co-expression in the same cells

I would add a limitations sentence such as:
> “Because all analyses were performed on bulk RNA-seq, we cannot determine whether ZP3 isoform usage arises from tumor cells, immune cells, or admixture effects.”

---

### B. FDR handling
#### Problem
You state BH correction in some places, but the reporting is not fully consistent about **what was corrected across what family**:
- tumor-vs-normal tests across transcripts
- mixed models across transcripts × features
- GSEA outputs
- multiple immune-feature correlations
- pan-cancer comparisons
- external cohort testing

For example, in the Results you frequently present raw P-values as if they were final inferential values, while elsewhere you say FDR was used. The reader cannot reconstruct the exact family of tests.

### Actionable fix
Define clearly in Methods:
- **what constitutes a family**
- whether FDR was across:
  - transcripts only
  - immune features only
  - all transcript-feature combinations
  - all cancer types
- whether the figures show **raw P** or **FDR q**

For tables and figures, report:
- effect size
- 95% CI where applicable
- raw P
- FDR q

For the main text, prefer **q-values** where multiple testing is involved.

---

### C. Cox proportional-hazards assumptions & ordinal grade encoding
#### Problem
You asked for Cox/grade handling, but the current manuscript **does not actually include a survival analysis** or Cox models. That is itself a concern: the title and framing suggest biomarker significance, but the study doesn’t test survival. If survival is omitted intentionally, fine—but then remove any survival-adjacent implications.

If you did run grade-encoding analyses elsewhere or intend to incorporate them, do **not** treat WHO grade as continuous/ordinal without justification. For glioma, grade is not a linear variable in a Cox model by default.

### Actionable fix
Either:
- **remove survival language entirely**, or
- add a proper survival analysis with:
  - Cox PH model
  - proportional hazards assumption checks
  - grade encoded as categorical unless a linear trend is explicitly justified
  - stratification by key molecular factors where appropriate

Given the current scope, I recommend **not adding survival claims** unless you can do them rigorously.

---

### D. Pan-cancer descriptive vs inferential claims
#### Problem
The manuscript uses “pan-cancer” language, but the 32 TCGA cancer-type analysis is still mostly **descriptive/meta-analytic** rather than inferential in a causal/generalizable sense. You cannot conclude that ZP3 is a pan-cancer immune marker from a pooled correlation with high heterogeneity (I²≈77%).

### Actionable fix
Use language like:
- “across-cancer descriptive consistency”
- “heterogeneous positive association”
- “broad but non-uniform pan-cancer pattern”

Avoid:
- “establishes”
- “validates”
- “generalizes across cancers” unless qualified by heterogeneity

---

### E. Threshold sensitivity in single-cell / short-read re-analysis
#### Problem
The external kallisto re-analysis is interesting, but the manuscript currently risks over-interpreting quantification instability as biology. With 51-bp single-end reads, isoform proportions are **highly threshold- and compatibility-dependent**. Your current explanation may be directionally correct, but you need to show that the result is not an artifact of:
- transcriptome reference choice
- effective-length correction
- kallisto settings
- transcript annotation version
- low-abundance filtering
- detection threshold definitions

### Actionable fix
Add sensitivity analyses for the external isoform re-analysis:
- alternative transcriptome annotation or release
- pseudoalignment with different parameters if reasonable
- show how FL/RI proportions change under detection thresholds
- report the number of samples in which both FL and RI exceed minimum abundance thresholds
- present a clear distinction between **“detected”** and **“quantified reliably”**

At minimum, soften the conclusion:
> “These data are consistent with read-length-dependent assignment bias”
rather than
> “systematic assignment bias”
unless you demonstrate it experimentally or with benchmarking.

---

## 3) Results reporting
### Major
**Section(s): Abstract; Results; Figure legends; Supplementary Tables**

### Problem 1: Numerical consistency
There are several places where the same analysis is described with slightly different numbers or emphasis.

Examples:
- Abstract: FL-high enriched for TNF-α/NF-κB signaling (NES=+2.45); Results says TNF-α signaling via NF-κB (same value, fine), but the number of pathways/depleted sets is summarized differently.
- Abstract: “six of seven immune features” after adjustment; Results and Methods agree, but the wording around “partial independence” should be tightened.
- Abstract: “sample-level correlations were modest, ρ=0.13–0.54”; Results says this is across GBM and LGG AP events, but the number of events tested is not always explicit.
- External cohort section: “all |ρ|<0.19, all P>0.4” is good, but later “minimum detectable |ρ| at α=0.05 is 0.35 and 0.40” should be tied to exact alpha and correction status.

### Actionable fix
Create a consistency pass:
- one table of **key headline numbers**
- one canonical version of each statistic used in Abstract, Results, and legends
- ensure **same N**, same tested comparison, same directionality, and same correction status everywhere

---

### Problem 2: Precision and overprecision
The manuscript often uses excessive precision for quantities that are not necessarily that stable:
- median proportions to 3 decimals
- ρ values to 2 decimals
- huge P-values to 1 significant digit

### Actionable fix
Use appropriate precision:
- medians: 2–3 decimals depending on scale
- correlations: 2 decimals
- P-values: scientific notation if small; otherwise exact to 2–3 sig figs
- confidence intervals: 2 decimals for correlations, not more

Do not imply false accuracy.

---

### Problem 3: Null/attenuated results need stronger, more even-handed handling
#### Problem
You do acknowledge the null external gene-level cohort result, which is commendable. But the Discussion still somewhat rescues it by emphasizing power limitations and by pivoting quickly to isoform re-analysis. That can read as if the null is only a nuisance obstacle.

### Actionable fix
Be more explicit that:
- the external **gene-level** replication was null
- this weakens claims of a broadly transferable ZP3–immune relationship
- the isoform-level internal evidence is suggestive, but not fully externally validated

I would recommend a sentence like:
> “The absence of external gene-level replication should temper any claim that ZP3 is a robust immune-context marker in glioma across cohorts.”

That will improve credibility.

---

## 4) Discussion restraint & attribution
### Major
**Section(s): Discussion; Introduction; Conclusion**

### Problem
The manuscript repeatedly approaches mechanistic attribution without direct evidence. In particular:
- “immune-context association”
- “myeloid-enriched immunosuppression”
- “TREM2-associated”
- “GPX4–ZP3 axis”
- “membrane-accessible ZP3”
- “cytoplasmic antigen”
- “systematically favours retained-intron transcript”

Some of these are literature-derived, but the manuscript should not merge prior mechanistic claims with your own association data unless clearly labeled as background hypothesis.

### Actionable fix
Add a firm distinction between:
1. **what is from the literature**
2. **what is inferred from your data**
3. **what is speculative**

For example:
- Do **not** imply your data establish a TREM2-linked mechanism.
- If the cited Cell 2026 GPX4–ZP3 study is mentioned, present it as **external background context**, not as validation of your dataset.
- Replace mechanistic phrasing with cautious language:
  - “consistent with”
  - “hypothesis-generating”
  - “may reflect”
  - “could be compatible with”

### Specific rewrite suggestion
The final paragraph of Discussion should not end with biological mechanistic next steps as if the association data already point to a defined pathway. Instead:
> “These data motivate targeted experiments to determine whether the observed FL/RI switch is cell-type specific and whether it has functional consequences in glioma immune biology.”

---

## 5) Figure–text–table consistency
### Major
**Section(s): Figures 1–5; Supplementary Tables S1–S4; Results**

### Problem
There are a few consistency gaps and one structural concern:

- Figure 4 legend states “Forest plot of mixed-effects coefficients … versus seven immune features,” but Supplementary Table S2 reports pooled sample-level correlations for only M2 and myeloid. The reader may assume all seven features were meta-analyzed.
- Figure 5 uses mapping rate and transcript proportions, but the text explains read-length dependence more strongly than the figure itself demonstrates.
- Supplementary Table S3 includes the external null plus diagnostics and LOCO/L2CO transportability references, which makes the table legend do too much.
- The title “isoform-resolved proxy for immune context” suggests a validated biomarker, but the actual figures show an association study with partial validation.

### Actionable fix
- Make sure every figure legend states exactly what is shown and what is not.
- If only M2 and myeloid underwent meta-analysis, say so in Results and Figure 4/Supp Table S2.
- Separate “external null gene-level assessment” from “internal robustness” in the supplementary material.
- Consider adding a small schematic figure summarizing the workflow and proxy boundary. That would make the manuscript easier to navigate.

---

## 6) Language, structure, and formatting
### Major
**Section(s): Title, Abstract, Introduction, Discussion, References**

### Problem
The manuscript is intelligible, but it is too long and self-referential for the amount of biology actually established. The structure sometimes reads like a methods-validation dossier rather than a focused biological paper.

Also, the reference list is extremely long for a focused report, and several references appear to be recent, tangential, or used mainly to support broad contextual statements. The bibliographic load is disproportionate to the core findings.

### Actionable fix
- Shorten the title.
- Reduce repetitive restatement of the “proxy boundary” across Abstract/Results/Discussion/Conclusion.
- Collapse some methodological detail into Supplementary Methods.
- Trim the Discussion to the core points:
  1. tumor-normal isoform shift
  2. modest sample-level proxy validity
  3. heterogeneous immune association
  4. platform dependence in short-read external re-analysis
  5. lack of external gene-level replication

### References
The reference list should be audited for:
- necessity
- relevance
- overcitation of very recent, possibly tangential glioma/microenvironment papers
- balance between classic methods and core biological claims

For a focused report, a very long reference list can signal compensatory citation rather than conceptual clarity. I would strongly recommend pruning background citations that do not directly support a statement in the manuscript.

---

## 7) Data/code availability & ethics statements
### Major
**Section(s): Data availability; Code availability; Ethics approval; AI-use disclosure**

### Problem
These statements are close to acceptable but not yet publication-grade.

### Actionable fix
You should include:
- a **persistent repository URL or DOI** now, not “upon acceptance,” if possible
- exact software versions and container/environment details
- whether any data are controlled-access and how they were handled
- whether the frozen result tables are public
- a statement that analyses are fully reproducible from public inputs

For ethics:
- “Ethical approval was therefore not required” is acceptable, but some journals prefer:
  > “This study used only public, de-identified data and did not require new ethics approval.”
- If any GEO datasets have their own consent/ethics conditions, note that they were used in accordance with those terms.

For AI-use:
- Good that it is disclosed.
- Ensure the statement complies with the target journal’s policy.
- Clarify that AI was used only for language editing, not scientific interpretation.

---

## 8) Overall narrative and contribution positioning vs the cited Cell 2026 GPX4–ZP3 study
### Major
**Section(s): Introduction; Discussion; Conclusion; References 15–16**

### Problem
The manuscript currently risks sounding like it is extending the Cell 2026 GPX4–ZP3 study mechanistically, but it actually does something different:
- that Cell study is mechanistic / pathway-oriented
- your study is transcript-usage / proxy-validation / cross-cohort association-oriented

This distinction must be made explicitly to avoid overclaiming novelty.

### Actionable fix
Position the contribution as follows:

**What this manuscript adds**
- transcript-usage resolution of ZP3
- comparison of FL vs retained-intron behavior
- cross-cancer association mapping
- explicit proxy-boundary framing
- external platform-dependent re-analysis showing why isoform-level transfer is difficult

**What it does not add**
- mechanistic proof of GPX4–ZP3 biology
- direct validation of ZP3 protein isoforms
- tumor microenvironment spatial localization
- functional evidence of immune modulation

### Suggested phrasing
> “In contrast to recent mechanistic work implicating GPX4–ZP3 biology, our study addresses a different question: whether transcript-abundance-derived ZP3 isoform proportions can serve as a bounded proxy for immune-context associations in bulk RNA-seq.”

That would sharpen novelty without overselling.

---

# Minor revision recommendations

## A. Abstract wording
### Minor
**Section: Abstract**

- “validated” appears too strong in places; replace with “evaluated,” “benchmarked,” or “tested.”
- “external isoform-level test” is good, but make clear it is a **platform-specific re-analysis**, not a universal validation.
- Consider shortening the Methods sentence; it is too dense.

---

## B. Terminology consistency
### Minor
**Section(s): throughout**

Use one consistent term for the same thing:
- either “full-length” or “FL”
- either “retained-intron” or “RI”
- either “proxy marker” or “proxy”
- either “ecological” or “cohort-level”

Avoid alternating between “PSI-like proportion,” “proportion measure,” and “isoform proportion” without defining which is the primary term.

---

## C. SpliceSeq terminology
### Minor
**Section(s): Results, Methods**

You say “alternative-promoter measurements” and “AP1 PSI” but the mapping between your FL isoform and the exact SpliceSeq event should be explained more explicitly. Which transcript structure corresponds to AP1 vs AP2? This needs a schematic or a brief explanatory note.

---

## D. Mixed-effects model details
### Minor
**Section: Methods**

Add:
- random intercept justification
- whether random slopes were considered
- whether residuals were checked
- whether predictors/outcomes were standardized
- whether cancer-type sample size imbalance affected estimation

---

## E. GSEA details
### Minor
**Section: Methods; Figure 3 legend**

Add the exact ranking metric and how ties were handled. Also state whether you used:
- Hallmark gene sets only
- size filtering thresholds
- permutation type

---

## F. External re-analysis wording
### Minor
**Section: Results; Figure 5; Discussion**

Avoid stating “most plausible explanation” too definitively unless backed by explicit benchmarking. Use:
- “consistent with”
- “suggests”
- “is compatible with”

---

## G. Supplementary tables
### Minor
**Section(s): Supplementary Tables S1–S4**

Add columns for:
- raw P
- FDR q
- effect size
- CI where possible
- exact N per test

---

## H. Reference list
### Minor
**Section: References**

The reference list should be trimmed and checked for relevance. Several 2026 references look tangential and may distract from the core manuscript. If the target journal permits, move some contextual citations to the Supplementary Discussion or reduce them substantially.

---

# Priority-ranked revision plan

## Highest priority
1. **Downgrade causal/biomarker language**; present the conclusion as a bounded association-level proxy.
2. **Clarify statistical families and FDR handling** across all tests.
3. **Tighten external validation claims**: internal robustness is not external validation.
4. **Resolve figure/text/table alignment** and numerical consistency.
5. **Add explicit bulk-data limitations** about cell-of-origin and mechanism.

## Medium priority
6. Shorten and simplify the title and Abstract.
7. Improve Methods transparency for mixed models, GSEA, and external kallisto re-analysis.
8. Strengthen discussion of null external gene-level replication.
9. Trim and focus the reference list.

## Lower priority
10. Language polishing and terminology standardization.
11. Add a workflow schematic if space permits.
12. Standardize precision and reporting format.

---

# Bottom line

This is a thoughtful and technically ambitious manuscript, but in its current form it **overstates the biological specificity and validation status** of the main finding. The most publishable version of the paper will be one that says:

- ZP3 FL proportion differs between tumor and normal tissue;
- it is ecologically concordant with an alternative-promoter signal;
- it shows modest sample-level agreement and heterogeneous pan-cancer association with immune features;
- it is sensitive to sequencing platform/read length;
- it is therefore a **bounded proxy**, not a validated mechanistic marker.

That framing would make the paper substantially stronger, more credible, and more likely to survive expert review.

If you want, I can next provide:
1. a **line-by-line “major edits to sentence-level wording”** version, or  
2. a **journal-style reviewer report** with “Overall score / Significance / Validity / Clarity / Recommendation.”