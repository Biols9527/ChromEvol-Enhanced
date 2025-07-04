# ChromEvol增强版：染色体进化分析技术文档

## 概述 (Overview)

此增强版祖先状态重建脚本融合了源于ChromEvol的复杂方法论。主要增强功能包括最大似然法、分支异质性建模、参数优化和随机映射。

**注意**：随着代码的更新（特别是 `src/ancestral_reconstruction.py` 中 `ChromEvolutionModel` 的修改），本文档中的某些细节（例如 `demidupl` 的具体实现和 `base_number` 的使用）可能与最新代码存在差异。请以代码实现和最新的 `README.md` 为准。本文档旨在提供一个关于 ChromEvol 核心思想如何被借鉴和扩展的概览。

## 主要特性 (Key Features)

### 1. 受ChromEvol启发的最大似然建模 (ChromEvol-Inspired Maximum Likelihood Modeling)

#### `ChromEvolutionModel` 类
- **参数化进化模型**，包含以下事件的速率：
  - `gain`: 单条染色体增加事件 (i -> i+1)
  - `loss`: 单条染色体丢失事件 (i -> i-1)
  - `dupl`: 全基因组复制事件 (i -> 2i)
  - `demidupl_gain`, `demidupl_loss`: 半复制后分别导致染色体数目增加或减少的事件。具体实现可能涉及多个目标状态。
  - `base_number_gain`, `base_number_loss`: 与基础数相关的染色体增减事件 (i -> i+b, i -> i-b)。
  - `linear_rate`: 染色体数目依赖的速率调整因子，影响`gain`和`loss`速率。
  - `base_number`: 模型中的基础染色体数目，用于`base_number_gain/loss`事件。

#### 转移速率矩阵 (Transition Rate Matrix)
- 遵循ChromEvol原理构建 **Q-矩阵** (瞬时速率矩阵)。
- 通过**矩阵指数运算**计算转移概率：P(t) = exp(Q*t)。
- 支持**分支特异性参数**以实现异质性模型。

### 2. 最大似然优化 (Maximum Likelihood Optimization)

#### `ChromEvolOptimizer` 类
- 使用 `scipy.optimize` 进行**参数优化**。
- 支持**多种优化算法** (如 L-BFGS-B)。
- **有界优化**确保速率参数为正。
- 提供**收敛诊断**和成功报告。

#### 模型选择 (Model Selection)
- 基于 **AIC/BIC标准**进行模型比较。
- 支持**多种模型类型**：例如，简单模型 vs. 完整模型。
- 根据信息准则**自动选择最佳模型**。

### 3. 基于似然的祖先状态重建 (Likelihood-Based Ancestral State Reconstruction)

#### `ChromEvolLikelihoodCalculator` 类
- 使用 **Felsenstein剪枝算法**计算似然性。
- **边际祖先状态重建**。
- 基于概率进行**最大似然状态分配**。
- 支持**根频率约束** (可指定根状态的先验信念)。

### 4. 随机映射 (Stochastic Mapping)

#### `StochasticMapper` 类
- 使用速率矩阵在分支上**模拟进化事件**。
- **蒙特卡洛采样** (默认：100次重复)。
- **事件类型分类**：增加、丢失、复制、基础数相关事件等。
- 计算期望事件数的**置信区间**。
- 生成**分支特异性的事件摘要**。

## 命令行界面 (Command Line Interface)

### 基本用法 (Basic Usage)

```bash
# 传统的最大简约法分析 (如果脚本支持此模式作为默认或通过参数切换)
python src/ancestral_reconstruction.py --tree <tree_file> --counts <counts_file> --out_image <image_file>

# 受ChromEvol启发的最大似然分析 (启用 --use_chromevol)
python src/ancestral_reconstruction.py --tree <tree_file> --counts <counts_file> --use_chromevol --out_image <image_file>

# 完整的优化与模型选择 (示例性，具体参数可能调整)
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --optimize_params \
    --model_selection \
    --out_image results/phylogeny_optimized.svg
```

### 主要参数 (Key Parameters)

#### ChromEvol特性参数
- `--use_chromevol`: 启用最大似然方法替代简约法。
- `--optimize_params`: 使用最大似然法优化模型参数 (需启用 `--use_chromevol`)。
- `--model_selection`: 执行基于AIC/BIC的模型选择 (需启用 `--use_chromevol`)。
- `--max_chr_number`: 模型中考虑的最大染色体数。
- `--root_count` (或类似参数如 `--root_prior_state`): 约束根节点的染色体数目 (在ML中用作先验)。

#### 输出参数 (增强)
- `--out_ancestors`: 输出的祖先状态报告，现在包含ML概率和期望事件数。
- `--out_image`: 增强的可视化输出，包含事件信息。

## 输出增强 (Output Enhancements)

### 祖先状态报告 (Ancestral States Report)
增强的CSV输出现在包含 (具体列名请参照最新脚本输出):
- `node_name`: 内部节点标识符。
- `inferred_chromosome_count`: 最可能的祖先状态。
- `ml_probability`: 该状态的最大似然概率。
- `expected_gain`, `expected_loss`, `expected_duplication`, 等: (来自随机映射的) 期望事件数。
- `state_probabilities`: (可选) 祖先节点各状态的完整后验概率分布。

### 增强的可视化 (Enhanced Visualization)
- 内部节点上显示**最大似然概率**。
- 分支上可能标注**期望事件数**。
- 可能用**橙色高亮**显示具有高复制活性的分支。
- **增强的图例**解释所有可视化元素。

### 控制台输出 (Console Output)
- **参数优化过程**，包括初始和最终参数值。
- 用于模型比较的**对数似然值**。
- **AIC/BIC模型选择**结果。
- **随机映射摘要**，包括各分支的期望事件数。

## 科学背景 (Scientific Background)

### ChromEvol的启发 (ChromEvol Inspiration)
此实现借鉴了ChromEvol方法论 (Mayrose et al. 2010, Glick & Mayrose 2014):

1.  **最大似然框架**: 使用概率模型而非简约法。
2.  **参数模型**: 为不同类型的染色体变化定义明确的速率。
3.  **分支异质性**: 允许不同分支使用不同的参数集。
4.  **随机映射**: 模拟进化历史以估计事件的期望数量。

### 数学框架 (Mathematical Framework)

#### 速率矩阵构建 (Rate Matrix Construction)
瞬时速率矩阵Q遵循ChromEvol的方法，但具体事件类型和参数化方式请参考 `ChromEvolutionModel` 类的最新实现。通常包括：
- `Q[i,i+1]` (染色体增加): 通常为 `gain_rate * (1 + linear_rate * i)`
- `Q[i,i-1]` (染色体丢失): 通常为 `loss_rate * (1 + linear_rate * i)`
- `Q[i,2*i]` (全基因组复制): `duplication_rate`
- 其他事件如半复制、基础数相关的增减，其速率贡献会加到相应的 `Q[i,j]` 上。
- 对角线元素 `Q[i,i] = -sum(Q[i,j] for j != i)`。

#### 似然计算 (Likelihood Calculation)
总似然性: L = Σ P(根状态=i) * L(数据|根状态=i)
其中 L(数据|根状态=i) 使用Felsenstein的剪枝算法计算。

#### 优化 (Optimization)
使用有界限的梯度下降等方法优化参数，以最大化对数似然函数。

## 性能考量 (Performance Considerations)

### 计算复杂度 (Computational Complexity)
- **矩阵指数运算**: O(n³)，其中n是最大染色体数。
- **随机映射**: O(重复次数 × 分支数 × 平均每分支事件数)。
- **优化**: 取决于收敛速度和参数空间的复杂性。

### 内存使用 (Memory Usage)
- 转移概率矩阵: 每个分支 O(n²)。
- 似然性数组: O(节点数 × n)。
- 随机映射: O(重复次数 × 分支数 × 事件数)。

## 验证与测试 (Validation and Testing)

### 测试案例 (Test Cases)
增强后的脚本已使用以下类型的树进行了测试：
- 小型树 (5-10 物种)
- 中型树 (20-50 物种)
- 大型树 (100+ 物种)
- 多种染色体数目范围 (例如, 5-50)

### 与ChromEvol的比较 (Comparison with ChromEvol)
与原始ChromEvol的主要区别可能包括：
- **模型细节**: 某些事件（如半复制）的模型定义可能简化或有所不同。
- **Python实现**: 相较于C++更易于访问和修改。
- **集成可视化**: 内建的树图绘制功能。
- **流线型接口**: 通常通过单个脚本执行。

## 未来增强方向 (Future Enhancements)

### 潜在扩展 (Potential Extensions)
1.  **基础数优化**: 自动推断基础染色体数目。
2.  **速率异质性**: 例如，跨分支的Gamma分布式速率。
3.  **模型复杂度**: 增加额外的转移类型 (如易位等)。
4.  **贝叶斯推断**: 使用MCMC采样评估参数的不确定性。
5.  **比较分析**: 多基因家族的染色体进化。

### 集成机会 (Integration Opportunities)
- **R集成**: 导出结果以便在R中进行详细统计分析。
- **可视化增强**: 使用Plotly等库实现交互式绘图。
- **数据库连接**: 直接访问染色体数据库。
- **流程集成**: 与Snakemake/Nextflow等工作流管理系统兼容。

## 参考文献 (References)

1.  Mayrose, I., Barker, M. S., & Otto, S. P. (2010). Probabilistic models of chromosome number evolution and the inference of polyploidy. *Systematic Biology*, 59(2), 132-144.
2.  Glick, L., & Mayrose, I. (2014). ChromEvol: assessing the pattern of chromosome number evolution and the inference of polyploidy along a phylogeny. *Molecular Biology and Evolution*, 31(7), 1914-1922.
3.  Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17(6), 368-376.

## 如何引用 (Citation)

如果您在研究中使用了此增强脚本，请引用：
- 原始ChromEvol论文 (Mayrose et al. 2010, Glick & Mayrose 2014)。
- 您对本增强实现的具体使用情况。

## 支持 (Support)

有关ChromEvol增强功能的疑问、bug报告或特性请求，请参考脚本的内联文档和注释，或在项目仓库中提交Issue。
