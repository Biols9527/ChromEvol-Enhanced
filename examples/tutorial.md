# 🎓 ChromEvol染色体进化分析教程

## 📚 教程概述

本教程将通过实际案例逐步指导您完成一个完整的染色体进化分析。我们将使用棘皮动物的数据集进行演示。

## 🎯 学习目标

完成本教程后，您将能够：
- ✅ 准备和验证输入数据
- ✅ 运行不同类型的祖先重建分析
- ✅ 解释分析结果
- ✅ 生成发表质量的可视化图形
- ✅ 进行比较分析和统计检验

## 📊 示例数据集介绍

### 研究背景
我们将分析后口动物(Deuterostomia)的染色体进化，特别关注：
- 棘皮动物各类群的染色体数目变化
- 祖先状态重建
- 染色体融合/分裂事件的定量分析

### 数据文件
- `pruned_phylogeny.nwk`: 包含39个物种的系统发育树
- `chromosome_counts.csv`: 各物种的染色体数目数据
- `all.map.tsv`: 共线性区块数据

## 🚀 步骤1: 环境准备和数据检查

### 1.1 检查工作目录
```bash
# 确认您在正确的目录中
pwd
ls -la

# 应该看到以下文件结构
# data/ docs/ src/ scripts/ results/ examples/
```

### 1.2 验证数据完整性
```bash
# 检查数据文件
ls data/
head -5 data/chromosome_counts.csv
head -3 data/pruned_phylogeny.nwk
```

**预期输出:**
```csv
species,chromosome_count
Antedon_bifida,11
Antedon_mediterranea,11
Leptometra_celtica,11
Branchiostoma_floridae,19
```

### 1.3 快速数据概览
```bash
# 统计物种数量
wc -l data/chromosome_counts.csv

# 查看染色体数目分布
cut -d',' -f2 data/chromosome_counts.csv | sort -n | uniq -c
```

## 🔬 步骤2: 基础分析 - 最大简约法

### 2.1 运行基础分析
```bash
# 使用默认的最大简约法
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --out_image results/tutorial_basic.svg \
    --out_ancestors results/tutorial_basic_ancestors.csv \
    --out_rearrangements results/tutorial_basic_events.csv
```

### 2.2 检查基础结果
```bash
# 查看祖先状态重建结果
head -10 results/tutorial_basic_ancestors.csv

# 查看重排事件
head -10 results/tutorial_basic_events.csv
```

### 2.3 结果解读
**关键观察点:**
- 根节点(Deuterostomia)的推断染色体数
- 主要类群的祖先状态
- 融合/分裂事件的分布

## 🧬 步骤3: 高级分析 - ChromEvol方法

### 3.1 运行ChromEvol分析
```bash
# 使用ChromEvol最大似然方法
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --map data/all.map.tsv \
    --use_chromevol \
    --out_image results/tutorial_chromevol.svg \
    --out_ancestors results/tutorial_chromevol_ancestors.csv \
    --out_rearrangements results/tutorial_chromevol_events.csv
```

### 3.2 参数优化分析
```bash
# 启用参数优化
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --optimize_params \
    --out_image results/tutorial_optimized.svg \
    --out_ancestors results/tutorial_optimized_ancestors.csv \
    --out_rearrangements results/tutorial_optimized_events.csv
```

### 3.3 模型选择分析
```bash
# 完整的模型选择和优化
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --optimize_params \
    --model_selection \
    --layout circular \
    --show_support \
    --out_image results/tutorial_complete.svg \
    --out_ancestors results/tutorial_complete_ancestors.csv \
    --out_rearrangements results/tutorial_complete_events.csv
```

## 📈 步骤4: 结果比较与解释

### 4.1 比较不同方法的结果
```bash
# 生成比较报告
python scripts/generate_summary.py

# 查看统计摘要
cat results/analysis_summary_table.csv
```

### 4.2 结果解释指南

#### 4.2.1 祖先状态置信度
```bash
# 查看ML概率分布
grep -v "node_name" results/tutorial_complete_ancestors.csv | \
cut -d',' -f3 | sort -n | tail -10
```

**解释标准:**
- ML概率 > 0.95: 高置信度，结果可靠
- ML概率 0.80-0.95: 中等置信度，需要谨慎解释
- ML概率 < 0.80: 低置信度，结果不确定

#### 4.2.2 进化事件分析
```bash
# 统计事件类型
grep "fusion" results/tutorial_complete_events.csv | wc -l
grep "fission" results/tutorial_complete_events.csv | wc -l
```

**生物学意义:**
- 融合事件多: 表明基因组趋向于紧缩
- 分裂事件多: 表明基因组趋向于扩张
- 事件大小: 反映进化压力强度

## 🎨 步骤5: 高级可视化和定制分析

### 5.1 不同布局比较
```bash
# 圆形布局 - 适合展示大型系统发育树
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --layout circular \
    --show_support \
    --show_branch_length \
    --out_image results/tutorial_circular.svg

# 矩形布局 - 适合详细分析
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --layout rectangular \
    --show_support \
    --show_branch_length \
    --out_image results/tutorial_rectangular.svg
```

### 5.2 物种对比分析
```bash
# 分析特定物种对之间的进化路径
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --analyze_pairs Antedon_bifida Leptometra_celtica Ophiocomina_nigra Asterias_rubens \
    --out_image results/tutorial_pairs.svg
```

### 5.3 高分辨率输出
```bash
# 生成发表质量的图片
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --optimize_params \
    --layout circular \
    --show_support \
    --dpi 300 \
    --out_image results/tutorial_publication.png
```

## 📊 步骤6: 数据分析和统计检验

### 6.1 运行综合统计分析
```bash
python scripts/generate_summary.py
```

### 6.2 关键统计指标解读

#### 模型拟合质量
```bash
# 查看模型选择结果
grep -A 10 "Model Selection" results/tutorial_complete_*.log
```

**评估标准:**
- AIC/BIC值越小越好
- 对数似然值越大(接近0)越好
- 参数收敛成功表明模型稳定

#### 进化速率估计
```bash
# 查看随机映射结果
grep "Expected Number of Events" -A 20 results/tutorial_complete_*.log
```

**解释要点:**
- 期望事件数反映分支的进化活跃度
- 数值越大表明该分支经历更多染色体重排
- 零值表明该分支相对稳定

## 🔍 步骤7: 结果验证和质量控制

### 7.1 生物学合理性检查
```python
# 创建验证脚本
cat > validate_results.py << 'EOF'
import pandas as pd

# 读取结果
ancestors = pd.read_csv('results/tutorial_complete_ancestors.csv')
events = pd.read_csv('results/tutorial_complete_events.csv')

# 检查染色体数目范围
chr_counts = ancestors['inferred_chromosome_count']
print(f"染色体数目范围: {chr_counts.min()} - {chr_counts.max()}")
print(f"平均染色体数目: {chr_counts.mean():.1f}")

# 检查置信度分布
ml_probs = ancestors['ml_probability'].dropna()
high_conf = (ml_probs > 0.95).sum()
med_conf = ((ml_probs > 0.80) & (ml_probs <= 0.95)).sum()
low_conf = (ml_probs <= 0.80).sum()

print(f"\n置信度分布:")
print(f"高置信度 (>0.95): {high_conf}")
print(f"中等置信度 (0.80-0.95): {med_conf}")
print(f"低置信度 (<0.80): {low_conf}")

# 检查事件类型分布
fusion_count = (events['type'] == 'fusion').sum()
fission_count = (events['type'] == 'fission').sum()
print(f"\n事件类型分布:")
print(f"融合事件: {fusion_count}")
print(f"分裂事件: {fission_count}")
print(f"融合/分裂比: {fusion_count/max(fission_count,1):.2f}")
EOF

python validate_results.py
```

### 7.2 与文献数据比较
- 比较推断的祖先状态与已知的化石记录
- 验证物种染色体数目与文献报道的一致性
- 检查进化趋势的生物学合理性

## 📝 步骤8: 报告撰写

### 8.1 结果摘要模板
```markdown
## 分析结果摘要

### 数据集
- 物种数量: 39
- 染色体数目范围: 11-25
- 系统发育覆盖: 后口动物主要类群

### 方法
- 祖先重建: ChromEvol最大似然法
- 模型选择: AIC/BIC准则
- 参数优化: 有界优化算法

### 主要发现
1. 后口动物祖先染色体数: 24条
2. 进化趋势: 融合占主导(比例4:1)
3. 海百合类群: 极端染色体缩减(45%减少)
4. 置信度: 主要节点高置信度(>95%)

### 生物学意义
- 染色体融合可能与基因组紧缩相关
- 不同类群表现出异质性的进化模式
- 海百合的特化可能反映了适应性进化
```

### 8.2 图表说明模板
```markdown
## 图1. 后口动物染色体进化系统发育树

圆形系统发育树显示了39个后口动物物种的染色体数目进化。
分支颜色表示进化事件类型：蓝色=融合，红色=分裂，绿色=稳定。
节点数字显示推断的祖先染色体数目和最大似然概率。
比例尺表示进化距离单位。
```

## 🎯 步骤9: 进阶分析选项

### 9.1 自定义参数分析
```bash
# 固定根节点染色体数
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --root_count 24 \
    --out_image results/tutorial_constrained.svg
```

### 9.2 批量分析脚本
```bash
# 创建批量分析脚本
cat > batch_analysis.sh << 'EOF'
#!/bin/bash

# 设置参数数组
layouts=("circular" "rectangular")
methods=("" "--use_chromevol" "--use_chromevol --optimize_params")

# 循环运行分析
for layout in "${layouts[@]}"; do
    for method in "${methods[@]}"; do
        echo "Running analysis: layout=$layout, method=$method"
        python src/ancestral_reconstruction.py \
            --tree data/pruned_phylogeny.nwk \
            --counts data/chromosome_counts.csv \
            --layout $layout \
            $method \
            --out_image results/batch_${layout}_${method// /_}.svg
    done
done
EOF

chmod +x batch_analysis.sh
./batch_analysis.sh
```

## 🏆 教程总结

### 完成的技能
- ✅ **数据准备**: 格式化和验证输入文件
- ✅ **基础分析**: 最大简约法祖先重建
- ✅ **高级建模**: ChromEvol最大似然方法
- ✅ **模型选择**: 统计模型比较
- ✅ **参数优化**: 进化参数估计
- ✅ **结果解释**: 生物学意义阐释
- ✅ **可视化**: 发表质量图形生成
- ✅ **质量控制**: 结果验证和检查

### 关键要点回顾
1. **方法选择**: 根据数据特点选择合适的分析方法
2. **模型验证**: 始终检查模型拟合质量和生物学合理性
3. **置信度评估**: 重视统计置信度，谨慎解释低置信结果
4. **多重验证**: 使用不同方法交叉验证结果
5. **文献比较**: 将结果与已知生物学知识对照

### 下一步建议
- 🔬 探索更大的数据集
- 📊 尝试其他类群的分析
- 🧬 结合基因组注释数据
- 📖 深入学习相关理论背景
- 🤝 与领域专家合作验证结果

---

🎉 **恭喜！您已经掌握了染色体进化分析的完整流程。现在可以应用这些技能到您自己的研究项目中了！**

如有疑问，请参考：
- 📚 [用户手册](user_guide.md)
- 🔬 [技术文档](ChromEvol_Enhancement_README.md)
- 📞 [获取帮助](../README.md#contact--support)
