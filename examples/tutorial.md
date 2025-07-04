# 🎓 ChromEvol增强版：染色体进化分析教程

## 📚 教程概述 (Tutorial Overview)

本教程将通过实际案例，一步步指导您完成一个完整的染色体进化分析流程。我们将以棘皮动物为例进行演示。

## 📊 示例数据集介绍 (Example Dataset Introduction)

### 研究背景 (Background)
我们将分析后口动物 (Deuterostomia) 的染色体进化，重点关注：
- 棘皮动物 (Echinodermata) 各主要类群的染色体数目变化。
- 重建关键祖先节点的染色体数目。
- 定量分析染色体融合 (fusion) 与分裂 (fission) 事件。

### 数据文件 (Data Files)
本项目 `data/` 目录下提供了示例数据：
- `pruned_phylogeny.nwk`: 包含39个物种的系统发育树文件（Newick格式）。
- `chromosome_counts.csv`: 包含各物种半倍体染色体数目的CSV文件。
- `all.map.tsv`: （可选）用于特定物种对之间染色体重排分析的共线性/同源性作图数据文件。

## 🚀 步骤1: 环境准备和数据检查 (Environment Setup and Data Check)

### 1.1 克隆仓库并安装依赖 (Clone Repository and Install Dependencies)
如果您尚未完成，请先克隆本项目的Git仓库，并安装所需的Python依赖包。
```bash
# 克隆仓库
# git clone https://github.com/Biols9527/ChromEvol-Enhanced.git
# cd ChromEvol-Enhanced

# 使用conda (推荐)
conda env create -n ChromEvol python==3.11 libgl
conda activate ChromEvol
pip install -r requirements.txt

```

### 1.2 检查工作目录和数据 (Verify Working Directory and Data)
确保您当前位于项目的主目录下。
```bash
pwd
ls -la data/
```
您应该能看到 `data/` 目录下的 `pruned_phylogeny.nwk` 和 `chromosome_counts.csv` 文件。

使用以下命令快速查看数据文件内容：
```bash
head -n 5 data/chromosome_counts.csv
head -n 3 data/pruned_phylogeny.nwk
```
**`chromosome_counts.csv` 预期输出示例:**
```csv
species,chromosome_count
Antedon_bifida,11
Antedon_mediterranea,11
Leptometra_celtica,11
Branchiostoma_floridae,19
```

### 1.3 数据概览 (Quick Data Overview)
对数据有一个初步的了解。
```bash
# 统计物种数量 (减去表头)
awk 'NR > 1' data/chromosome_counts.csv | wc -l

# 查看染色体数目分布的概况
awk -F',' 'NR > 1 {print $2}' data/chromosome_counts.csv | sort -n | uniq -c
```

## 🔬 步骤2: 基础分析 - 最大简约法 (Basic Analysis - Maximum Parsimony)

最大简约法 (Maximum Parsimony, MP) 是一种较为简单的祖先状态重建方法，它寻找需要最少进化改变次数的进化路径。

### 2.1 运行基础简约法分析
我们的脚本 `ancestral_reconstruction.py` 在不指定 `--use_chromevol` 时，会采用简约法（或一种简约法思想的实现）进行重建。
```bash
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --out_image results/tutorial_parsimony.svg \
    --out_ancestors results/tutorial_parsimony_ancestors.csv \
    --out_rearrangements results/tutorial_parsimony_events.csv \
    --root_count 24 \
    --max_chr_number 60 # 假设最大染色体数不超过60
```
- `--root_count 24`: 我们根据文献假设后口动物的共同祖先有24条染色体。这对简约法是一个强约束。
- `--max_chr_number 60`: 设置模型能处理的最大染色体数。

### 2.2 检查简约法结果
查看生成的CSV文件：
```bash
head -n 10 results/tutorial_parsimony_ancestors.csv
head -n 10 results/tutorial_parsimony_events.csv
```
同时，您可以在 `results/` 目录下找到 `tutorial_parsimony.svg` 图像文件，用图像浏览器打开查看。

### 2.3 简约法结果解读
- **祖先状态**: `tutorial_parsimony_ancestors.csv` 中列出了内部节点的推断染色体数。
- **进化事件**: `tutorial_parsimony_events.csv` 中记录了各分支上发生的融合或分裂事件。
- **关键观察点**:
    - 根节点 (Deuterostomia) 的推断染色体数是否符合预期 (24)。
    - 主要进化枝（如海星、海胆、海参等）的祖先状态。
    - 融合与分裂事件在系统发育树上的分布模式。

简约法的结果提供了一个基础的演化图景，但它没有考虑不同事件发生的概率或分支长度。

## 🧬 步骤3: 高级分析 - ChromEvol启发的最大似然法 (Advanced Analysis - ChromEvol-Inspired Maximum Likelihood)

最大似然法 (Maximum Likelihood, ML) 基于概率模型来推断祖先状态和进化参数。`--use_chromevol` 参数会启用此方法。

### 3.1 运行基本的ChromEvol ML分析 (不优化参数)
使用模型默认参数进行分析。
```bash
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --root_count 24 \
    --max_chr_number 60 \
    --out_image results/tutorial_chromevol_basic.svg \
    --out_ancestors results/tutorial_chromevol_basic_ancestors.csv \
    --out_rearrangements results/tutorial_chromevol_basic_events.csv
```
- `--use_chromevol`: 启用ML方法。
- 输出文件增加了 `_chromevol_basic` 后缀以区分。

### 3.2 启用参数优化 (Optimize Parameters)
允许脚本优化模型的速率参数（如 gain, loss, duplication 等速率）以更好地拟合数据。
```bash
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --optimize_params \
    --root_count 24 \
    --max_chr_number 60 \
    --out_image results/tutorial_chromevol_optimized.svg \
    --out_ancestors results/tutorial_chromevol_optimized_ancestors.csv \
    --out_rearrangements results/tutorial_chromevol_optimized_events.csv
```
- `--optimize_params`: 启用参数优化。
- 分析过程会稍慢，因为需要进行迭代优化。控制台会输出优化过程信息。

### 3.3 进行模型选择 (Model Selection)
比较不同复杂度的模型（例如，只包含gain/loss的模型 vs 包含所有事件类型的完整模型），并根据AIC/BIC信息准则选择最优模型。
```bash
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --optimize_params \
    --model_selection \
    --root_count 24 \
    --max_chr_number 60 \
    --layout circular \
    --show_branch_length \
    --show_support \
    --out_image results/tutorial_chromevol_modelselection.svg \
    --out_ancestors results/tutorial_chromevol_modelselection_ancestors.csv \
    --out_rearrangements results/tutorial_chromevol_modelselection_events.csv
```
- `--model_selection`: 启用模型选择。脚本会测试多个预定义的模型配置。
- `--layout circular`: 将输出图像设置为圆形布局。
- `--show_branch_length` 和 `--show_support`: 在图中显示分支长度和支持率（如果树文件包含这些信息）。

运行此命令后，注意观察控制台输出，会显示不同模型的AIC/BIC值以及最终选择的模型。后续的祖先状态重建和随机映射将基于所选的最佳模型进行。

## 📈 步骤4: 结果解读与比较 (Interpreting and Comparing Results)

### 4.1 生成综合摘要报告
使用 `scripts/generate_summary.py` 脚本来汇总分析结果。
```bash
python scripts/generate_summary.py \
    --ancestral_states results/tutorial_chromevol_modelselection_ancestors.csv \
    --rearrangement_events results/tutorial_chromevol_modelselection_events.csv \
    --chromosome_counts data/chromosome_counts.csv \
    --output_summary_table results/tutorial_summary_table.csv
```
- 根据您上一步选择分析的输出文件名，相应修改 `--ancestral_states` 和 `--rearrangement_events` 参数。

该脚本会在控制台打印详细的摘要信息，并将一个简洁的摘要表保存到 `results/tutorial_summary_table.csv`。

### 4.2 结果解读指南

#### 4.2.1 祖先状态与置信度
查看 `_ancestors.csv` 文件 (例如 `tutorial_chromevol_modelselection_ancestors.csv`)：
- `inferred_chromosome_count`: 推断的祖先染色体数目。
- `ml_probability`: 对该推断状态的最大似然概率（或后验概率，取决于具体实现）。值越高，置信度越高。
- `expected_...` 列: 随机映射得到的各类期望事件数。

**评估标准示例:**
- ML概率 > 0.95: 高置信度，结果较为可靠。
- ML概率 0.80 - 0.95: 中等置信度，解释时需谨慎。
- ML概率 < 0.80: 低置信度，结果不确定性较大。

#### 4.2.2 进化事件分析 (Evolutionary Event Analysis)
查看 `_events.csv` 文件 (例如 `tutorial_chromevol_modelselection_events.csv`)，此文件基于推断的祖先状态直接计算染色体数目变化：
- `event_type`: 'fusion' (融合) 或 'fission' (分裂)。
- `change_magnitude`: 变化的染色体数目。

查看 `generate_summary.py` 的控制台输出或其生成的摘要表：
- **融合/分裂事件总数和比例**: 了解哪种类型的事件在进化中占主导。
- **随机映射期望事件数**: 例如 `总期望gain事件数`, `总期望loss事件数`, `总期望duplication事件数` 等，这些提供了更细致的事件类型概率估计。

**生物学意义联想:**
- **融合事件为主**: 可能表明基因组趋向于紧缩，或某些特定的选择压力。
- **分裂事件为主**: 可能与基因组扩张、新基因产生或适应性辐射有关。
- **高频率的全基因组复制 (duplication) 事件**: 可能指示了多倍化的关键节点，对物种形成和多样性有重要影响。

## 🎨 步骤5: 高级可视化与定制分析 (Advanced Visualization and Custom Analysis)

### 5.1 不同布局与显示选项
尝试不同的可视化参数组合：
```bash
# 矩形布局，显示分支长度和支持率
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol --optimize_params \
    --root_count 24 --max_chr_number 60 \
    --layout rectangular \
    --show_branch_length \
    --show_support \
    --out_image results/tutorial_rectangular_detailed.svg
```

### 5.2 特定物种对之间的重排分析 (Synteny-based Rearrangement Analysis)
如果提供了 `--map` 文件 (如 `data/all.map.tsv`)，可以分析指定物种对之间的染色体结构变异。
```bash
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --map data/all.map.tsv \
    --use_chromevol \
    --analyze_pairs Antedon_bifida Leptometra_celtica Ophiocomina_nigra Asterias_rubens \
    --out_image results/tutorial_pairs_analysis.svg
    # 其他参数如 --root_count, --max_chr_number 等可以按需添加
```
- `--analyze_pairs`: 后跟一系列物种名，脚本会两两配对进行分析。
- 控制台会输出这些物种对之间的融合/分裂事件（基于map文件推断）。

### 5.3 生成高分辨率图像用于出版
使用 `--dpi` 参数（主要对PNG格式有效）。对于SVG或PDF这类矢量格式，它们本身是无损缩放的。
```bash
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol --optimize_params --model_selection \
    --root_count 24 --max_chr_number 60 \
    --layout circular \
    --dpi 300 \
    --out_image results/tutorial_publication.png
```

## 📊 步骤6: 深入解读统计指标 (Interpreting Key Statistical Metrics)

回顾 `scripts/generate_summary.py` 生成的控制台输出和CSV表格。

#### 模型拟合与参数 (Model Fit and Parameters)
- **对数似然值 (Log-likelihood)**: 在参数优化和模型选择过程中会输出。值越大（越接近0），模型对数据的拟合越好。
- **AIC/BIC值**: 用于模型选择。值越小，模型相对越优（在拟合优度和参数数量之间取得平衡）。
- **优化的参数值**: 例如 `gain_rate`, `loss_rate` 等。这些值反映了不同类型染色体数目变化事件的相对速率。例如，如果 `duplication_rate` 显著高于其他速率，则表明全基因组复制可能是该类群进化的一个重要模式。

#### 随机映射的期望事件数 (Expected Event Counts from Stochastic Mapping)
- `ancestral_states_report.csv` 中各内部节点（分支）的 `expected_...` 列，以及 `generate_summary.py` 输出的总期望事件数。
- 这些数值反映了在系统发育树的特定分支上，平均预期会发生多少次特定类型的事件。
- “进化最活跃的谱系”部分指出了哪些分支经历了更多的染色体数量变化。

## 🔍 步骤7: 结果验证与质量控制 (Results Validation and Quality Control)

### 7.1 生物学合理性检查 (Biological Plausibility Check)
- **推断的祖先染色体数是否合理？** 例如，一个深层祖先节点的染色体数不太可能远超其所有后代。
- **主要进化事件发生的位置是否符合已有的生物学认知？** 例如，如果已知某个类群经历过近期多倍化，随机映射是否检测到了相应的复制信号？
- **染色体数目变化的趋势是否与类群的生活史、基因组大小等特征相关联？**

### 7.2 与文献数据比较 (Comparison with Literature)
- 将推断的关键祖先状态与文献中基于其他证据（如细胞学、比较基因组学）的研究进行比较。
- 验证本分析中使用的物种染色体数目与权威数据库或文献报道是否一致。

### 7.3 检查日志输出
仔细阅读分析过程中控制台（或日志文件，如果配置了）的 `WARNING` 和 `ERROR` 信息，这些可能指示数据问题或分析中的潜在陷阱。

## 📝 步骤8: 撰写分析报告 (Composing the Analysis Report)

一个典型的分析报告可能包含以下部分：

### 8.1 数据集与方法摘要
- **数据集**: 物种数量，染色体数目范围，系统发育树的来源和覆盖范围。
- **方法**:
    - 祖先重建方法 (简约法 vs 最大似然法)。
    - 如果使用ML：所采用的进化模型（例如，参数化的ChromEvol模型，包含哪些事件类型），是否进行了参数优化和模型选择，随机映射的设置。

### 8.2 主要发现 (Key Findings)
1.  **关键祖先节点的染色体数目**: 例如，目标类群的共同祖先 (MRCA) 的推断染色体数及其置信度/概率。
2.  **主要的进化趋势**: 染色体数目是倾向于增加（分裂为主）还是减少（融合为主）？是否存在显著的全基因组复制事件？
3.  **特定进化枝的模式**: 某些子类群是否表现出与其他类群不同的进化模式（例如，某个分支经历了快速的染色体数目变化）？
4.  **模型参数和速率**: 优化得到的各事件类型的相对速率如何？

### 8.3 生物学意义探讨 (Biological Interpretation)
- 观察到的染色体进化模式可能与哪些生物学过程相关？（例如，物种形成速率、适应性辐射、基因组大小变化、生活史策略的转变等）。
- 分析结果是否支持或挑战了关于该类群进化的现有假说？

### 8.4 图表与说明
- **图1: 注释的染色体进化系统发育树**。图注应清晰说明各可视化元素的含义（颜色编码、节点标签等）。
- **表1: 关键祖先节点染色体数目及支持度**。
- **表2: (可选) 模型比较结果 (AIC/BIC值)**。
- **表3: (可选) 优化后的模型参数**。

## 🎯 步骤9: 进阶分析选项与探索 (Advanced Options and Exploration)

### 9.1 自定义模型参数进行分析
如果您对某些进化速率有先验知识，或者想测试特定假设，可以在不进行优化的前提下，手动指定模型参数来运行分析。这通常需要在代码层面修改 `ChromEvolutionModel` 的默认参数或通过特定参数传入（如果脚本支持）。

### 9.2 批量分析与比较不同约束的影响
使用脚本自动化运行多种参数组合（例如，不同的根节点约束、不同的模型复杂度），然后比较结果的稳健性。

```bash
# 示例 batch_analysis.sh
# #!/bin/bash
#
# # 定义要测试的根节点计数值
# root_counts=(20 24 28)
#
# for rc in "${root_counts[@]}"; do
#     echo "Running analysis with root_count = $rc"
#     python src/ancestral_reconstruction.py \
#         --tree data/pruned_phylogeny.nwk \
#         --counts data/chromosome_counts.csv \
#         --use_chromevol \
#         --optimize_params \
#         --root_count $rc \
#         --max_chr_number 60 \
#         --out_image results/tutorial_root_${rc}.svg \
#         --out_ancestors results/tutorial_root_${rc}_ancestors.csv \
#         --out_rearrangements results/tutorial_root_${rc}_events.csv
# done
```
您可以比较不同根节点约束下的结果差异。

## 🏆 教程总结 (Tutorial Summary)

### 关键要点回顾 (Key Takeaways)
1.  **方法选择**: 根据研究问题和数据特性，选择合适的分析方法（简约法 vs 最大似然法）。ML方法通常更灵活且信息更丰富，但计算量也更大。
2.  **模型验证与比较**: 如果使用ML方法，进行模型选择（如AIC/BIC比较）非常重要，以避免模型过于简单或过度参数化。
3.  **参数的生物学意义**: 理解模型中各参数（如gain, loss, duplication速率）的含义，并结合优化结果进行解释。
4.  **结果的稳健性**: 测试不同参数（特别是如根节点状态这样的强假设）对结果的影响。
5.  **整合多方证据**: 将本流程的分析结果与文献报道、其他类型的基因组数据（如比较基因组学、细胞遗传学）结合，进行综合性的生物学解释。

### 下一步建议 (Next Steps)
- 🔬 在您自己的研究数据上应用本流程。
- 📊 深入探索不同模型配置（例如，在 `ChromEvolutionModel` 中尝试不同的事件组合或约束）。
- 🧬 如果您的研究涉及特定基因家族的进化，思考如何将染色体进化背景与基因层面的进化关联起来。
- 📖 进一步阅读染色体进化、多倍化以及相关生物信息学方法的文献。

---
