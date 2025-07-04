# 📖 用户使用手册

## 🎯 概述

本用户手册将指导您完成染色体进化分析的全部流程，从数据准备到结果解释。

## 🛠️ 安装指南

### 步骤1: 下载项目
```bash
git clone https://github.com/Biols9527/ChromEvol-Enhanced.git
cd ChromEvol-Enhanced
```

### 步骤2: 安装依赖
```bash
# 使用conda (推荐)
conda env create -n ChromEvol python==3.11 libgl
conda activate ChromEvol
pip install -r requirements.txt
```

### 步骤3: 验证安装
```bash
export export QT_QPA_PLATFORM=offscreen
or
QT_QPA_PLATFORM=offscreen python src/ancestral_reconstruction.py --help
```

## 📁 数据准备

### 必需的输入文件

#### 1. 系统发育树文件 (Newick格式)
```
# 示例: pruned_phylogeny.nwk
((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05);
```

#### 2. 染色体数目文件 (CSV格式)
```csv
species,chromosome_count
Species_A,22
Species_B,24
Species_C,20
Species_D,18
```

#### 3. 共线性数据文件 (TSV格式) - 可选
```tsv
species_A	chromosome_A	species_B	chromosome_B	mapping_type_A_to_B	mapping_type_B_to_A	bidirectional_mapping_type	count
Amphiura_filiformis	chr5	Ophiocomina_nigra	OZ221938.1	1:1	1:1	1:1	932
Amphiura_filiformis	chr3	Ophiocomina_nigra	OZ221941.1	1:2	1:1	1:n (Fission)	679
```

### 数据质量检查
- 确保物种名称在所有文件中一致
- 检查染色体数目的合理性 (通常在1-100之间)
- 验证系统发育树的格式正确性

## 🚀 基本使用

### 快速开始
```bash
# 最简单的分析
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv
```

### 标准分析
```bash
# 使用ChromEvol方法的完整分析
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --map data/all.map.tsv \
    --use_chromevol \
    --optimize_params \
    --model_selection \
    --layout circular \
    --show_support \
    --out_image results/analysis.svg
```

## 🔧 高级参数设置

### 分析方法选择
```bash
# 使用最大简约法 (默认)
python src/ancestral_reconstruction.py --tree tree.nwk --counts counts.csv

# 使用ChromEvol最大似然法
python src/ancestral_reconstruction.py --tree tree.nwk --counts counts.csv --use_chromevol

# 启用参数优化
python src/ancestral_reconstruction.py --use_chromevol --optimize_params

# 进行模型选择
python src/ancestral_reconstruction.py --use_chromevol --model_selection
```

### 可视化设置
```bash
# 圆形布局
--layout circular

# 矩形布局  
--layout rectangular

# 显示分支长度
--show_branch_length

# 显示支持值
--show_support

# 设置图片分辨率
--dpi 300
```

### 输出控制
```bash
# 指定输出文件
--out_image results/my_tree.svg
--out_ancestors results/ancestral_states.csv
--out_rearrangements results/events.csv

# 约束根节点染色体数
--root_count 24
```

### 物种对比分析
```bash
# 分析特定物种对之间的重排事件
--analyze_pairs Species_A Species_B Species_C Species_D
```

## 📊 结果解读

### 输出文件说明

#### 1. 祖先状态重建结果 (*_ancestral_states.csv)
```csv
node_name,inferred_chromosome_count,ml_probability,expected_gains,expected_losses,expected_duplications
Root,24,0.95,0.1,0.05,0.02
Node1,22,0.88,0.15,0.1,0.01
```

**字段含义:**
- `node_name`: 节点名称
- `inferred_chromosome_count`: 推断的染色体数目
- `ml_probability`: 最大似然概率
- `expected_gains/losses/duplications`: 期望的事件数

#### 2. 重排事件报告 (*_rearrangements.csv)
```csv
branch,type,change,from_count,to_count
Branch_A,fusion,5,24,19
Branch_B,fission,2,20,22
```

**字段含义:**
- `branch`: 分支名称
- `type`: 事件类型 (fusion/fission)
- `change`: 变化量
- `from_count/to_count`: 起始/结束染色体数

#### 3. 系统发育树图片
- **蓝色分支**: 融合事件 (染色体数减少)
- **红色分支**: 分裂事件 (染色体数增加)
- **绿色分支**: 稳定 (染色体数不变)
- **橙色分支**: 高倍增活动
- **灰色虚线**: 缺失或不确定数据

### 统计指标解释

#### 模型选择标准
- **AIC (赤池信息准则)**: 数值越小，模型越好
- **BIC (贝叶斯信息准则)**: 数值越小，模型越好
- **对数似然值**: 数值越大 (越接近0)，拟合越好

#### 置信度评估
- **ML概率 > 0.95**: 高置信度
- **ML概率 0.80-0.95**: 中等置信度  
- **ML概率 < 0.80**: 低置信度

## 🔍 故障排除

### 常见错误及解决方案

#### 1. 导入错误
```
ImportError: No module named 'ete3'
```
**解决方案**: 重新安装依赖包
```bash
pip install ete3 pandas numpy scipy
```

#### 2. 数据格式错误
```
ValueError: Species names mismatch between tree and counts file
```
**解决方案**: 检查并统一所有文件中的物种名称

#### 3. 内存不足
```
MemoryError: Unable to allocate array
```
**解决方案**: 
- 减少物种数量
- 增加系统内存
- 使用更简单的模型

#### 4. 优化失败
```
OptimizationError: Parameter optimization failed to converge
```
**解决方案**:
- 使用 `--optimize_params` 时尝试不同的初始参数
- 检查数据质量
- 使用简约法作为备选方案

### 性能优化建议

#### 大数据集处理
```bash
# 对于大型数据集，推荐使用简约法
python src/ancestral_reconstruction.py --tree large_tree.nwk --counts large_counts.csv

# 关闭参数优化以加快速度
python src/ancestral_reconstruction.py --use_chromevol  # 不使用 --optimize_params
```

#### 内存优化
- 使用SVG格式输出 (相比PNG更节省内存)
- 关闭不必要的可视化选项
- 分批处理大型数据集

## 📝 最佳实践

### 数据预处理
1. **质量控制**: 验证染色体数目的生物学合理性
2. **物种名标准化**: 使用一致的命名规范
3. **系统发育树校验**: 确保拓扑结构正确

### 分析流程
1. **探索性分析**: 先用默认参数运行
2. **模型选择**: 比较不同模型的拟合效果
3. **参数优化**: 在确定模型后进行参数优化
4. **结果验证**: 检查生物学合理性

### 结果报告
1. **记录参数**: 保存完整的命令行参数
2. **统计检验**: 报告置信度和显著性
3. **生物学解释**: 结合生物学背景解释结果
