# -*- coding: utf-8 -*-
"""Apply bridge-dual-review v0.5 revision replacements to Article3_A3English_draft_v0.5.md."""
import io

PATH = r"C:\D\workbuddy\Research\Extracellular_GPX4_Immunosuppression\article3\manuscripts\Article3_A3English_draft_v0.5.md"

pairs = [
# R1 Results Fig2 proxy marker -> proxy definition
('''We therefore describe the FL proportion as a proxy marker of ZP3 splicing that is valid at the cohort/ecological scale but not a validated substitute for sample-level splice-event measurement.''',
'''We therefore describe the FL proportion as a proxy of ZP3 splicing—an association-level, non-causal measure—that is concordant at the cohort/ecological scale but not a validated substitute for sample-level splice-event measurement.'''),

# R2 Results Fig4 section: add bulk inference boundary
('''### The FL proportion associates with immune features across 32 cancer types (Fig. 4)

In mixed-effects models [25] across 9,186 tumor samples''',
'''### The FL proportion associates with immune features across 32 cancer types (Fig. 4)

All associations in this section are sample-level and derived from bulk RNA-seq; they do not identify the expressing cell type and do not establish co-expression within the same cells.

In mixed-effects models [25] across 9,186 tumor samples'''),

# R3 meta-analysis: descriptive cross-cancer pattern caveat
('''To move beyond the GBM/LGG-only sample-level evidence, we performed a cancer-stratified fixed-effect meta-analysis [75] of Spearman(FL, immune score) within each of the 32 TCGA cancer types with ≥30 samples.''',
'''To move beyond the GBM/LGG-only sample-level evidence, we performed a cancer-stratified fixed-effect meta-analysis [75] of Spearman(FL, immune score) within each of the 32 TCGA cancer types with ≥30 samples. We emphasize that this constitutes a descriptive cross-cancer pattern with substantial heterogeneity, not evidence of uniform generalization across cancers.'''),

# R4 external gene-level null: strengthen transparent framing
('''We report this null result transparently.''',
'''We report this null result transparently: the external gene-level replication of the ZP3–immune association was null, which weakens any claim that ZP3 is a robust immune-context marker across cohorts and defines a boundary on generalization.'''),

# R5 external isoform: systematically favours -> read-length-dependent assignment effects
('''transcript-compatibility assignment therefore systematically favours the retained-intron transcript, which contains all exonic and intronic sequence.''',
'''transcript-compatibility assignment therefore favours the retained-intron transcript under this read-length regime, which contains all exonic and intronic sequence; we interpret this as a read-length-dependent assignment effect rather than evidence of systematic quantification bias, although alternative annotations or thresholds could shift the proportions (see Methods).'''),

# R6 external isoform closing: platform-specific + detected vs reliably quantified
('''This analysis adds a seventh, external step to our evaluation workflow and sharpens the proxy boundary: transcript-abundance-derived isoform proxies should be validated within the same read-length/sequencing platform on which they will be applied.''',
'''This analysis adds a seventh, external step to our evaluation workflow and sharpens the proxy boundary: transcript-abundance-derived isoform proxies should be evaluated within the same read-length/sequencing platform on which they will be applied. We further distinguish "detected" (ZP3 total TPM > 0, 24/24 samples) from "reliably quantified" at isoform resolution, which was constrained by the 51-bp single-end design.'''),

# R7 Discussion: define proxy as non-causal
('''This distinction matters because unvalidated proxy measures can produce reproducible but misleading associations [14,53,54].''',
'''This distinction matters because unvalidated proxy measures can produce reproducible but misleading associations [14,53,54]. We use "proxy" as a non-causal, context-associated measure: it carries no implication of cell-of-origin, protein-level validation, or mechanism.'''),

# R8 Discussion: most plausible -> most parsimonious + caveat
('''The most plausible explanation is mechanical rather than biological: short single-end reads frequently fail to span the junctions that define the full-length transcript, and pseudo-alignment assigns such reads to the transcript with the greatest compatible sequence—here, the retained-intron transcript, which includes the entire intronic sequence.''',
'''The most parsimonious interpretation is mechanical rather than biological—short single-end reads frequently fail to span the junctions that define the full-length transcript, and pseudo-alignment assigns such reads to the transcript with the greatest compatible sequence—although we cannot exclude contributions from annotation choice, effective-length correction, or low-abundance quantification noise (see Methods).'''),

# R9 Discussion limitations: bulk inference boundary
('''Limitations include the use of transcript-level proportions rather than junction-level PSI [60,61]; modest sample-level agreement with independent splice measurements; substantial cross-cancer heterogeneity in the sample-level signal; the absence of protein-level validation of isoform usage; and the observational design [62–67].''',
'''Limitations include the use of transcript-level proportions rather than junction-level PSI [60,61]; modest sample-level agreement with independent splice measurements; substantial cross-cancer heterogeneity in the sample-level signal; the absence of protein-level validation of isoform usage; and the observational design [62–67]. Because all analyses were performed on bulk RNA-seq, we cannot determine whether ZP3 isoform usage arises from tumor cells, immune cells, or admixture effects, and the associations do not establish co-expression within the same cells.'''),

# R10 Conclusion: claim calibration
('''Because sample-level agreement with independent splice-event measurements is modest, we present the FL proportion as an explicitly bounded proxy marker rather than a validated isoform replacement, and we offer the multi-scale evaluation workflow as a reusable template for transcript-abundance-derived isoform proxies.''',
'''Because sample-level agreement with independent splice-event measurements is modest, the external gene-level association was null, and isoform-level composition is read-length dependent, we present the FL proportion as an explicitly bounded, association-level proxy—not a validated isoform replacement and not a mechanistic marker—and we offer the multi-scale evaluation workflow as a reusable template for transcript-abundance-derived isoform proxies.'''),

# R11 Methods: expand statistical reporting into inference & multiple testing
('''### Statistical reporting and reproducibility

Tests were two-sided unless stated otherwise; FDR used Benjamini–Hochberg; P values were not interpreted as evidence of causality.''',
'''### Statistical inference and multiple testing

All units of analysis are samples (patients or tissues), not cells; immune-feature scores are bulk RNA-seq-derived signatures/deconvolution outputs, and all isoform–immune associations are sample-level. Tests were two-sided unless stated otherwise; FDR used the Benjamini–Hochberg procedure. Multiple-testing families were defined per analysis module: (i) tumor-versus-normal transcript comparisons (corrected across the seven transcripts); (ii) mixed-effects immune-feature associations (corrected across transcripts × immune features); (iii) GSEA (FDR q from the pre-ranked algorithm); (iv) compositional-control and meta-analysis correlations (nominal P, reported with effect sizes and 95% CIs); and (v) external cohort tests (nominal P per cohort, pooled 95% CI). Main-text associations that involve multiple testing are reported with q-values or with P plus explicit correction status; where nominal P is shown (external cohorts, compositional controls), it is labelled as such and should not be read as confirmatory. Exact N varies with filtering (e.g., 7,577 after low-signal filtering) and is stated at each analysis. P values were not interpreted as evidence of causality.'''),

# R12 Methods kallisto: detected vs reliably quantified caveat
('''All 24 samples passed pseudo-alignment (representative mapping rate 87.2%; SRR7050184: 56.8 M reads, 49.5 M pseudo-aligned).''',
'''All 24 samples passed pseudo-alignment (representative mapping rate 87.2%; SRR7050184: 56.8 M reads, 49.5 M pseudo-aligned). Because isoform proportions from short single-end reads are sensitive to transcriptome annotation, effective-length correction, and detection thresholds, we report "detected" (ZP3 total TPM > 0) separately from "reliably quantified" at isoform resolution, and we interpret the external proportions as platform-specific estimates rather than definitive measurements.'''),

# R13 Figure 4 legend: meta scope clarification
('''Forest plot of mixed-effects coefficients (cancer type as random effect; n=9,186 tumor samples) for FL (top, positive) and retained-intron (bottom, negative) proportions versus seven immune features. FL: M2 β=0.28 (P=8.4×10−35), myeloid β=0.22, IFN-γ β=0.20, Treg β=0.12, checkpoint β=0.12, T-cell exhaustion β=0.11, cytolytic n.s. RI: reciprocal negative. The cross-cancer pooled sample-level correlation (meta-analysis) is reported in Supplementary Table S2. Source: `a3_mixed_model_frozen.csv`.''',
'''Forest plot of mixed-effects coefficients (cancer type as random effect; n=9,186 tumor samples) for FL (top, positive) and retained-intron (bottom, negative) proportions versus seven immune features. All associations are sample-level and derived from bulk RNA-seq; they do not identify the expressing cell type. FL: M2 β=0.28 (P=8.4×10−35), myeloid β=0.22, IFN-γ β=0.20, Treg β=0.12, checkpoint β=0.12, T-cell exhaustion β=0.11, cytolytic n.s. RI: reciprocal negative. The cancer-stratified meta-analysis of the pooled sample-level correlation was performed for M2 macrophages and myeloid cells only (Supplementary Table S2); the remaining features were analysed in the mixed-effects models shown here. Source: `a3_mixed_model_frozen.csv`, `a3_robustness_meta.csv`.'''),

# R14a Figure 5 legend title: platform-specific
('''### Figure 5. External isoform-level re-analysis of GSE113474 (n=24 GBM) with kallisto.''',
'''### Figure 5. Platform-specific isoform-level re-analysis of GSE113474 (n=24 GBM) with kallisto.'''),

# R14b Figure 5 legend closing: exploratory caveat
('''(c) Pseudo-alignment mapping rate by sample (representative 87.2%). Source: `a3_external_zp3_isoform.csv`, `a3_external_kallisto_tpm.csv`, `a3_external_validation_summary.md`.''',
'''(c) Pseudo-alignment mapping rate by sample (representative 87.2%). This re-analysis demonstrates read-length/platform dependence of isoform assignment; it is an exploratory, platform-specific assessment, not a full quantification benchmark across annotations. Source: `a3_external_zp3_isoform.csv`, `a3_external_kallisto_tpm.csv`, `a3_external_validation_summary.md`.'''),

# R15 Data availability: repository wording
('''A public, versioned copy with a persistent identifier (DOI) will be deposited in a public repository (e.g., Zenodo) upon acceptance; the exact URL and commit hash will be inserted at submission.''',
'''A public, versioned copy with a persistent identifier (DOI) will be deposited in a public repository (e.g., Zenodo) prior to publication; the exact URL and commit hash will be provided at submission and updated upon deposit.'''),

# R16 Ethics: standard wording
('''This study used only publicly available, de-identified datasets and did not involve new human-participant recruitment, animal experiments or access to identifiable clinical records. Ethical approval was therefore not required.''',
'''This study used only public, de-identified datasets and did not require new ethics approval; it involved no new human-participant recruitment, no animal experiments, and no access to identifiable clinical records. GEO datasets were used in accordance with their terms and consent conditions.'''),
]

with io.open(PATH, "r", encoding="utf-8") as f:
    text = f.read()

ok, miss, multi = [], [], []
for i, (old, new) in enumerate(pairs, 1):
    n = text.count(old)
    if n == 0:
        miss.append(i)
    elif n == 1:
        text = text.replace(old, new, 1)
        ok.append(i)
    else:
        multi.append((i, n))

with io.open(PATH, "w", encoding="utf-8", newline="\n") as f:
    f.write(text)

print(f"OK: {len(ok)} items -> {ok}")
print(f"MISSING: {len(miss)} -> {miss}")
print(f"MULTI(should be none): {multi}")
print(f"total replacements attempted: {len(pairs)}")
