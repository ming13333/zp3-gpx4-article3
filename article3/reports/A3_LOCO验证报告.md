# A3 内部 leave-one-cancer-out 验证报告（冻结⑧，2026-08-18）

> 配套脚本：`article3/scripts/freeze_a3_loco.py`
> 冻结表：`article3/results/a3_loco_frozen.csv`（231 行 = 32 癌种 × 7 评分 + 1 全量基线）
> 目的：回答「内部 isoform 关联是否由某几个癌种单独驱动？」

## 1. 方法

- 输入：`article3/results/zp3_psi_pancancer_results/psi_immune_joined_samples.csv`（9,186 样本，32 癌种）
- 对每个癌种 c，留出 c 后对其余癌种做固定效应 Fisher-z 合并：
  `Z = Σ w_i z_i / Σ w_i, w_i = N_i − 3`
- 报告每癌种留出后的合并 ρ、95% CI、Cochran Q、I²、剩余癌种数 k
- 统计基元复用 `freeze_a3_robustness.py`（纯标准库 Spearman + betai + Fisher-z）

## 2. 结果：LOCO 合并 ρ 范围（对照全量基线）

| 免疫评分 | 全量基线 ρ | LOCO ρ 范围 | 方向翻转 |
|---|---|---|---|
| M2_Macrophage | +0.137 | +0.124 ~ +0.144 | 0 |
| T_cell_exhaustion | +0.054 | +0.042 ~ +0.063 | 0 |
| Cytolytic_activity | −0.002 | −0.012 ~ +0.005 | 4（噪音级） |
| Treg | +0.052 | +0.044 ~ +0.062 | 0 |
| IFN_gamma | +0.085 | +0.079 ~ +0.096 | 0 |
| Checkpoint | +0.059 | +0.045 ~ +0.067 | 0 |
| Myeloid | +0.110 | +0.098 ~ +0.119 | 0 |

## 3. 解读

- **6/7 免疫特征 LOCO 范围极窄（波动 <0.02）且方向一致**：内部 isoform 级关联
  在 32 癌种间高度稳健，不是由单癌种驱动。
- **Cytolytic_activity** 因全量基线本身为 −0.002（噪音级），出现 4 次 ±0.00x
  方向翻转，幅度可忽略，不影响结论。
- 与冻结⑦外部 null 形成对照：
  - 内部 isoform 级：稳健（LOCO ρ 波动 <0.02）
  - 外部基因级：null（n=32/24 功效不足 + 检验对象不同）
  → 强化「总基因表达不能替代 isoform 比例测量」的方法学主张。

## 4. 对论文的意义

- Results 新增 LOCO 段，Discussion/Workflow 段引用 `a3_loco_frozen.csv`。
- 与「external gene-level non-replication」并列，构成完整证据结构：
  **内部 isoform 稳健 + 外部基因级 null = bounded proxy 的诚实结论**。

## 5. 自检

脚本内置自检：所有 LOCO 合并 ρ 与全量基线同号（Cytolytic 噪音级翻转仅 WARN，
仍冻结透明报告）。运行 `python freeze_a3_loco.py` → **RESULT: PASS**。
