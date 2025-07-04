# 🧬 ChromEvol增强版：染色体进化分析流程

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](docs/)

> **一个受ChromEvol方法启发的、用于染色体进化分析的先进生物信息学流程**

## 🎯 项目概述 (Overview)

本项目提供了一个全面的工具包，用于分析系统发育树中的染色体进化，主要功能包括最大似然法祖先状态重建、染色体重排分析以及生成出版质量的可视化结果。该流程集成了源于ChromEvol的复杂方法与现代Python科学计算技术。

## ✨ 主要特性 (Key Features)

- 🔬 **受ChromEvol启发的最大似然建模**
- 📊 **祖先染色体状态重建**
- 🧬 **染色体重排量化** (融合、分裂、复制事件)
- 📈 **随机映射** 用于进化事件估计
- 🎨 **出版级的系统发育树可视化**
- 🔧 **灵活的命令行界面**
- 📋 **全面的统计报告**

## 🚀 快速开始 (Quick Start)

### 安装 (Installation)

```bash
# 克隆仓库
git clone https://github.com/Biols9527/ChromEvol-Enhanced.git
cd ChromEvol-Enhanced

# 安装依赖
pip install -r requirements.txt
```

### 基本用法 (Basic Usage)

```bash
# 使用ChromEvol方法运行完整分析
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --map data/all.map.tsv \
    --use_chromevol \
    --optimize_params \
    --model_selection \
    --out_image results/phylogeny.svg

# 生成综合摘要报告
python scripts/generate_summary.py \
    --ancestral_states results/ancestral_states_report.csv \
    --rearrangement_events results/rearrangement_events_report.csv \
    --chromosome_counts data/chromosome_counts.csv \
    --output_summary_table results/analysis_summary_table.csv
```

## 📁 项目结构 (Project Structure)

```
chromevol-enhanced-analysis/
├── 📂 src/                              # 源代码
│   └── ancestral_reconstruction.py      # 主分析流程脚本
├── 📂 scripts/                          # 工具脚本
│   └── generate_summary.py              # 统计摘要生成脚本
├── 📂 data/                             # 输入数据文件
│   ├── pruned_phylogeny.nwk            # 系统发育树 (Newick格式)
│   ├── chromosome_counts.csv           # 物种染色体数目数据
│   └── all.map.tsv                     # 同源性/共线性作图数据
├── 📂 results/                          # 分析输出结果
│   ├── *.csv                           # 统计报告
│   └── *.svg, *.png, *.pdf             # 可视化图片
├── 📂 docs/                             # 文档
│   ├── README.md                       # 技术文档 (ChromEvol增强功能)
│   ├── user_guide.md                   # 用户指南 (待补充)
│   └── api_reference.md                # API参考 (待补充)
├── 📂 examples/                         # 示例分析
│   └── tutorial.md                     # 分步教程
├── 📄 README.md                        # 本文件 (项目总览)
├── 📄 requirements.txt                 # Python依赖包
└── 📄 LICENSE                          # 许可证信息
```

## 🔧 环境要求 (Requirements)

- Python 3.7+
- 所需的Python包 (会自动安装):
  - `ete3` - 系统发育树处理
  - `pandas` - 数据分析
  - `numpy` - 数值计算
  - `scipy` - 科学计算与统计建模
  - `matplotlib` - 基础绘图

## 📖 文档 (Documentation)

- 📚 **[用户指南](docs/user_guide.md)** - 详细的使用说明 (待补充)
- 🔬 **[技术文档](docs/ChromEvol_Enhancement_README.md)** - 实现细节 (部分内容待翻译)
- 📊 **[科学报告示例](docs/evolutionary_analysis_report.md)** - 分析方法和结果示例 (待补充)
- 🎓 **[教程](examples/tutorial.md)** - 分步操作示例 (中文)

## 🎨 输出示例 (Example Output)

该流程可以生成出版质量的可视化结果，具有以下特点:
- 分支颜色编码 (蓝色=融合, 红色=分裂, 绿色=稳定)
- 统计注释 (最大似然概率, 支持率等)
- 多种布局选项 (圆形/矩形)
- 高分辨率输出格式 (SVG, PNG, PDF)

## 📊 支持的分析类型 (Supported Analysis Types)

### 1. 祖先状态重建
- 最大似然估计
- 最大简约法推断
- 统计置信度评估

### 2. 染色体重排分析
- 融合/分裂事件检测
- 事件量级量化
- 系统发育趋势分析

### 3. 模型选择与优化
- AIC/BIC 模型比较
- 参数优化
- 速率异质性建模

### 4. 随机映射
- 事件速率估计
- 蒙特卡洛模拟
- 分支特异性统计

## 🤝 贡献代码 (Contributing)

我们欢迎各种贡献！请查看我们的 [贡献指南](CONTRIBUTING.md) (英文) 获取详细信息。

1. Fork 本仓库
2. 创建您的特性分支 (`git checkout -b feature/amazing-feature`)
3. 提交您的更改 (`git commit -m 'Add amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 开启一个 Pull Request

## 📄 许可证 (License)

本项目采用MIT许可证授权 - 详情请见 [LICENSE](LICENSE) 文件。

## 🙏 致谢 (Acknowledgments)

- [ChromEvol](https://github.com/anatshafir1/chromevol) 的开发者们，感谢其方法论的启发
- ETE3 开发团队，感谢其提供的系统发育树处理工具
- Scientific Python 社区，感谢其提供的基础库支持

---

**⭐ 如果您觉得本项目有用，请考虑给它一个星标！**
