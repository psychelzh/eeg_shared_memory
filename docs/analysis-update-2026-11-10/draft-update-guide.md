# 完成分析后的 draft 更新清单

本次使用 2026 年 11 月 10 日 14:58 完成的 HPC 缓存，重新执行 `index.qmd` 的后续统计、模型拟合和绘图。主报告为 `_site/index.html`，当前报告的 16 张结果图已重新生成至 `figures/`，冻结输出同步更新至 `_freeze/index/`。本清单核对的是 `draft/individualized_memory_eeg_20260513.docx`；该文件中已包含主要补充图图注。独立的 `individualized_memory_eeg_supplementary.docx` 仅有少量占位文字，不是当前完整的补充结果稿。

下面的 P 编号为该 Word 文件正文 XML 中的段落序号，包含空段落。实际修改时用小节标题及引文关键词定位更方便。Word 原稿未直接改写。

## 最先处理的结论与表述

1. **移除仍残留的中介或因果结论。** Abstract P13、Introduction P38、Discussion P132、P139、P142 仍有 “mediated”、通过提高 fidelity 促进记忆等措辞。当前 Figure 3 是联合回归与共同性分解，支持 ISS 和 IGS 的共享及独特解释量，不能据此声称中介机制得到检验。P94–98 的模型数值与当前结果基本一致，可保留。
2. **区分全时段和滑动窗结果。** 同时控制语义与记忆能力后，SMC 全时段三个成分的 Holm 校正结果均不显著；但 C1 在 566–723 ms 的动态效应保留。P125、Abstract 和 Discussion 应明确这一层次，避免写成全部分析均支持独立效应。
3. **群体表征控制后未发现正向显著 SMC 效应。** 当前控制 group representation 后的分析不支持“超越群体表征仍可预测 SMC”的泛化结论。`draft/文章的主要结果.docx` P33 等旧笔记也应同步核对。
4. **不要把样本内拟合称为样本外预测。** Figure 3A 的预测值来自在同一批 206 人上拟合的联合模型，并非交叉验证。建议使用 “model-fitted memory scores” 或明确 “in-sample”。

Abstract / Discussion 可采用的关联性表述：

> Semantic alignment and group-aligned representational fidelity explained both shared and distinct portions of individual differences in memory performance. Pairwise representational similarity remained associated with shared memory content in a late C1 time window after controlling for semantic representations and memory ability.

## Behavioral results 和 Figure 1 与 Figure S1

定位：P46、P48、Figure 1 图注 P50，以及 Figure S1 图注 P237。

| 指标 | 当前结果 |
|---|---:|
| Overall d′，均值 ± SD | 1.333 ± 0.567 |
| Recollection d′，均值 ± SD | 1.355 ± 0.572 |
| Familiarity d′，均值 ± SD | 0.791 ± 0.530 |
| Overall 与 recollection | r = 0.845 |
| Overall 与 familiarity | r = 0.702 |
| Recollection 与 familiarity | r = 0.416 |
| SMC，均值 ± SD | 0.653 ± 0.056 |

- P46 的 “r’s > 0.78” 应删除，改报上表的实际相关；特别是 familiarity 与 overall 的相关为 0.702。
- Recollection 的 skew / kurtosis 更新为 −0.017 / −0.054；familiarity 为 0.593 / 0.573。Overall 的描述基本不变。
- P48 和 P50 的 memorability–individual content 平均相关由 0.278 更新为 **0.332**。若“平均解释方差”按每位被试相关的平方再取平均定义，则为 **12.33%**，原稿的 7.7% 应更新。`mean(r^2)` 与 `mean(r)^2` 不相同，不要混用。
- Figure 1 的任务示意图为手工组图部分；替换统计面板时保留任务示意图和既有设计信息。

来源：[行为描述](memory_descriptives.csv)、[行为相关](memory_correlations.csv)、[SMC 与 memorability](additional_behavior_descriptives.csv)。

## ISS 与记忆以及 Figure 2

定位：“Greater semantic representations predict better memory” P60、P62、P65。

- P60 的 ISS 单样本检验更新为 C1/C2/C3：**t(205) = 9.407 / 7.019 / 13.642**，均为 Holm p < .001。旧值 7.95、6.33、12.1 不再适用。
- 全时段 ISS–overall memory 的相关维持 C1/C2/C3：r = 0.159 / 0.237 / 0.327，Holm p = .0228 / .00123 / .00000475，可按原稿精度报告。
- ISS 成分比较中，C3–C1 的调整后 p 为 **.000861**，原稿 .008 应更正。这里的成分两两比较来自 `emmeans`，当前代码默认使用 **Tukey** 调整，不能标记为 Holm。C1–C2 的调整后 p 为 .0447，C3–C2 < .001。相关及聚类图中的 Holm 调整与这一比较方法要分开说明。
- P62 的 C1 508–664 ms 仅为未通过 Holm 的候选簇，当前 p_perm = .01399、p_Holm = .06993；C3 215–332 ms 的 p_perm = .02997、p_Holm = .11988。正文的显著结论应聚焦 **C3 410–684 ms，p_Holm = .005994**。
- P65 还保留了“逐人回归、R² 后取 log”的旧聚合方法。当前 Figure 3 使用 **三个成分的 ISS 求和得到 cISS、IGS 求和得到 cIGS**，与 Methods P181 一致，应改写并消除 “Figure SXX” 占位。

## IGS 与记忆以及 Figure 2

定位：P83、P85、P87。

- 全时段 IGS 的单样本检验及与 overall memory 的相关与原稿基本一致；C1、C3 相关显著，C2 不显著。
- P87 的动态时间窗必须更新：

| 成分 | 当前显著时间窗 | p_perm | p_Holm |
|---|---:|---:|---:|
| C1 | 566.4–800.8 ms | .000999 | .008991 |
| C3 | 449.2–878.9 ms | .000999 | .008991 |

原稿 C1 470–626 ms、C3 373–763 ms 不再适用。Methods 和 Discussion 中引用这些时间窗的句子也应同步改。

## 稳定性 GSS 和 Figure S2 与 Figure S4

定位：P74、P76、P78，Figure S2 图注 P241，Figure S4 图注 P246。

- 实际 split-half 抽样为**每组 20–100 人的两组互不重叠被试**。每组 100 人时，两组总计 200 人；当前平均稳定性 C1/C2/C3 为 **0.532 / 0.412 / 0.467**。原稿 “三成分均超过 0.5” 及 0.68 / 0.69 / 0.55 的数值不能作为该样本规模下的实测值报告。
- Figure S2 上排现仅在实际观测的每组 20–100 人范围显示拟合曲线，已移除向每组 200 人的外推；下排 GSS 仍使用实际单组 20–200 人的结果。两种横轴含义应在图注中分清。
- P76 的全时段 GSS 保持 C1/C2/C3：r = .139 / .136 / .401，三者 Holm p = .0003，原稿的总体结论可保留。
- P78 的 GSS 动态结果替换如下；在当前补充图中，GSS 动态位于 **Figure S4F**，原文指向 S4B 应更正。

| 成分 | 当前显著时间窗 ms | p_perm | p_Holm |
|---|---:|---:|---:|
| C1 | 293.0–449.2 | .006993 | .034965 |
| C1 | 546.9–722.7 | .001998 | .013986 |
| C2 | 156.3–312.5 | .004995 | .029970 |
| C2 | 605.5–898.4 | .000999 | .008991 |
| C3 | 214.8–644.5 | .000999 | .008991 |

## Figure 3 的独特和共享贡献

定位：图注 P91，Results P94–98，Methods P180–181，以及 Abstract / Discussion 的相应解释。

这部分数值已重算，原稿数值基本正确：cISS–cIGS r = .461；联合模型 R² = .15495，F(2, 203) = 18.611；cISS 标准化 β = .274，95% CI [.130, .417]；cIGS β = .184，95% CI [.040, .327]，p = .0122。

| 方差部分 | R² | 占总解释量 |
|---|---:|---:|
| ISS 独特 | .0589 | 38.04% |
| IGS 独特 | .0266 | 17.17% |
| 共享 | .0694 | 44.79% |

保留 P94–98 的主要数值；把残留中介结论统一改为联合回归与共同性分解所支持的共享及独特关联。Figure 3A 的拟合值与观测值相关约 .394，是样本内结果。旧中文结果笔记中的 .04 / .02 / .10 是另一聚合方案的数值，不应与当前 sum 方案混合。

来源：[联合模型](joint_model.csv)、[回归系数](joint_model_coefficients.csv)、[聚合分数相关](joint_score_correlations.csv)、[共同性分析](commonality.csv)。

## ISPS 和 SMC 以及 Figure 4 与 Figure S6

定位：P119、P121、P123、P125、P128、P143 和补充图注 P254。

- P119 的 ISPS 均值维持 C1/C2/C3：.008525 / .004755 / .004278，但三者当前 **p_Holm = .002997**，不能再写 p_Holm < .001。
- P121 的 C2 第二个 ISPS 时间窗改为 **585.9–761.7 ms**，原稿约 605–762 ms 不准确。其余显著窗及约 645 / 234 / 156 ms 的峰值基本保持；完整峰值见 [temporal_peaks.csv](temporal_peaks.csv)。
- P123 全时段 ISPS–SMC 相关为 .03267 / .02364 / .03039，当前对应 Holm p 为 **.0375 / .0398 / .0375**。
- C1 的 **566.4–722.7 ms** 动态 SMC 簇保留，当前 p_perm = .000999，**p_Holm = .006993**，原稿 .002 / .014 应更新。P123 的时间序列图引用应为 Figure 4E，Figure 4C 是 ISPS 本身的动态结果。
- 控制记忆能力后同一 C1 动态窗保留，p_Holm = .006993。
- 控制语义后同一动态窗保留，p_Holm = .005994；**同时控制语义与记忆能力**后仍保留，p_Holm = .006993，见 Figure S6F。
- 同时控制语义与记忆能力后的**全时段**相关为 C1/C2/C3：r = .02989 / .02218 / .01943，Holm p = **.0621 / .1026 / .1026**，均未显著。P125 应具体写动态 C1，而非不加区分地说 residual representations 预测 SMC。
- 控制 group representation（以及同时控制 memory ability）后没有通过校正的正向动态簇；Figure S6C–D 不支持“独立于群体表征”的泛化结论。

可用于 P125 的替换句：

> After controlling for semantic representations and individual memory ability, residual intersubject representational similarity remained associated with shared memory content in C1 during 566–723 ms (cluster-permutation p = .000999, Holm-adjusted p = .006993). The corresponding whole-epoch associations did not survive Holm correction.

## Recollection familiarity 与区域补充分析

定位：P60、P62、P87 的泛化表述，Figure S3 图注 P244，以及新增区域图。

- 全时段：C3 的 ISS 与 IGS 对 recollection、familiarity 均仍显著；C1/C2 并非全部显著。不要写成所有成分都重复主分析。
- ISS–recollection 动态显著窗：C2 **605.5–898.4 ms**、C3 **429.7–683.6 ms**，两者 p_Holm = .01998。
- ISS–familiarity：显著窗为 **C1 390.6–664.1 ms，p_Holm = .006993**；C3 在当前 Holm 校正后没有显著簇。因此旧稿 “similar temporal patterns” 应收窄。
- IGS–recollection：C1 **566.4–722.7 ms，p_Holm = .03996**；C3 动态簇均未通过 Holm。
- IGS–familiarity：C1 **312.5–488.3 ms，p_Holm = .03996** 和 **566.4–781.3 ms，p_Holm = .00999**；C3 **468.8–585.9 ms，p_Holm = .03596**。
- 三张六脑区图已更新，应在补充材料中分配正式编号，并说明 Region 1–6 对应 `config/eeg_channel_labels_64.csv` 中的通道定义。当前 overall、recollection、familiarity 的区域记忆相关及时间簇均未通过图中使用的 Holm 校正；这些结果不能写成显著区域定位证据。

## Methods 和图注需要同步的规则

定位：Behavioral scoring、P181、Mantel P205、cluster-permutation P207，以及 Figure 2 / Figure S4 / Figure S6 图注。

- Overall d′ 使用二分类 log-linear 校正 `(count + 0.5)/(N + 1)`。Remember/Know 使用四类别每类 0.5 的平滑，即分母 N + 2；familiarity 再按 `know/(1-rem)` 调整后计算 signal–noise 的 d′。原稿公式需要与当前 `R/behav.R` 一致。
- cISS / cIGS 使用三个成分的原始分数求和，回归时再标准化。不要保留 P65 中另一套 R²-log 聚合方法。
- Mantel 保留 9,999 次内部置换；时间簇使用 1,000 次外部置换。置换矩阵时，同一重排应作用于行列标签，不能描述为仅打乱一边行而不对应处理列。
- 置换 p 值采用加一校正。外层 p 为 `(b+1)/(1000+1)`，最小约 .000999；内部 Mantel p 最小为 .0001。经过 Holm 后的值应如实报告，不能沿用旧种子产生的 p。
- P207 的 “maximum cluster statistics over all six regions” 与当前流程不符：当前按 CCA 成分或区域分别形成时间簇及对应零分布，再在各输出 target 的簇 p 上进行 Holm 调整。图注应说明实际校正范围。
- P207 不能把所有 cluster mass 都称为 Fisher’s Z 之和：SMC–Mantel 动态使用 Mantel r 求和，ISPS 动态使用 ISPS 均值求和，其他检验按各自的统计量计算。
- `emmeans` 的 ISS/IGS 成分两两比较为 Tukey 调整；散点图的三成分相关及时间簇图使用 Holm。不要把全部分析统一标成 Holm。

## 图片替换对应

| 稿件位置 | 新结果文件 |
|---|---|
| Figure 1 统计面板 | `figures/fig-1-behav.png` |
| Figure 2 统计面板 | `figures/fig-2-iss-igs-pred-mem.png` |
| Figure 3 | `figures/fig-3-iss-igs-unique-shared.png` |
| Figure 4 统计面板 | `figures/fig-4-isps-and-smc.png` |
| Figure S1 | `figures/fig-supp-1-behav-corr.png` |
| Figure S2B–C | `figures/fig-supp-2-group-averaged-stability.png` |
| Figure S3 | `figures/fig-supp-3-iss-igs-pred-mem.png` |
| Figure S4 | `figures/fig-supp-4-iss-igs-comparison.png` |
| Figure S5 | `figures/fig-supp-5-word-shape-pattern.png` |
| Figure S6 | `figures/fig-supp-6-isps-smc-partial.png` |
| 反应锁定控制 | `figures/fig-supp-X-iss-igs-pred-mem-rt-locked.png` |
| 残差表征指标 | `figures/fig-supp-X-iss-igs-unique-partial.png` |
| 不同聚合方法的共同性分析 | `figures/fig-supp-X-iss-igs-commonality-robust.png` |
| 六脑区 overall | `figures/fig-supp-X-iss-mem-regions.png` |
| 六脑区 recollection | `figures/fig-supp-X-iss-igs-regions-recollection.png` |
| 六脑区 familiarity | `figures/fig-supp-X-iss-igs-regions-familiarity.png` |

`X` 图号需要在稿件中正式编号。历史的 `fig-3-iss-igs-mem.png`、`fig-X-iss-igs-unique-shared.png` 和 mediation 图片不是本次报告引用的新版图，不要误用。任务示意图、CCA 静态头皮地形图以及 Illustrator 手工组图需在替换统计面板时保留或另行组装。

## 可复查数据与复现

- [全时段统计](whole_time_statistics.csv)：68 个当前 target，含未校正 p 与该 target 内的 Holm p。
- [全部时间簇](temporal_clusters.csv) / [Holm 显著簇](significant_temporal_clusters.csv)：79 个当前 cluster target，包含精确毫秒时间、原始 p、Holm p 和置换次数。空的 start/end 表示未形成候选簇，不是缺失计算。
- [成分比较](component_comparisons.csv) / [split-half 稳定性](group_stability.csv)。
- 其余行为、模型、共同性分析 CSV 位于本目录；所有表均来自当前完成的缓存，没有重新启动昂贵置换。
- 重新导出：在项目根目录运行 `Rscript scripts/summarize_completed_results.R`。
- 重新计算报告并作图：`quarto render index.qmd --to html --execute`。文档级显式执行用于刷新冻结缓存。

本清单不是对论文全部方法及结论的完整审计；未单独验证原稿中 CCA 成分回归系数差值及 bootstrap CI 等没有直接出现在这次导出中的陈述。此类数值应追溯对应分析记录后再保留，不能用联合模型的 β 或 CI 替代。
