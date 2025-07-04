# 🧬 ChromEvol-Enhanced Chromosome Evolution Analysis

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](docs/)

> **A state-of-the-art bioinformatics pipeline for chromosome evolution analysis inspired by ChromEvol methodology**

## 🎯 Overview

This project provides a comprehensive toolkit for analyzing chromosome evolution across phylogenetic trees, featuring maximum likelihood ancestral state reconstruction, chromosomal rearrangement analysis, and publication-quality visualizations. The pipeline integrates sophisticated methodologies inspired by ChromEvol with modern Python scientific computing.

## ✨ Key Features

- 🔬 **ChromEvol-inspired maximum likelihood modeling**
- 📊 **Ancestral chromosome state reconstruction**
- 🧬 **Chromosomal rearrangement quantification** (fusions, fissions, duplications)
- 📈 **Stochastic mapping** for evolutionary event estimation
- 🎨 **Publication-ready phylogenetic visualizations**
- 🔧 **Flexible command-line interface**
- 📋 **Comprehensive statistical reporting**

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-username/chromevol-enhanced-analysis.git
cd chromevol-enhanced-analysis

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```bash
# Run complete analysis with ChromEvol methods
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --map data/all.map.tsv \
    --use_chromevol \
    --optimize_params \
    --model_selection \
    --out_image results/phylogeny.svg

# Generate comprehensive summary
python scripts/generate_summary.py
```

## 📁 Project Structure

```
chromevol-enhanced-analysis/
├── 📂 src/                              # Source code
│   └── ancestral_reconstruction.py      # Main analysis pipeline
├── 📂 scripts/                          # Utility scripts
│   └── generate_summary.py              # Statistical summary generator
├── 📂 data/                             # Input data files
│   ├── pruned_phylogeny.nwk            # Phylogenetic tree
│   ├── chromosome_counts.csv           # Species chromosome data
│   └── all.map.tsv                     # Synteny/collinearity data
├── 📂 results/                          # Analysis outputs
│   ├── *.csv                           # Statistical reports
│   └── *.svg, *.png, *.pdf             # Visualizations
├── 📂 docs/                             # Documentation
│   ├── README.md                       # Technical documentation
│   ├── user_guide.md                   # User manual
│   └── api_reference.md                # API documentation
├── 📂 examples/                         # Example analyses
│   └── tutorial.md                     # Step-by-step tutorial
├── 📄 README.md                        # This file
├── 📄 requirements.txt                 # Python dependencies
└── 📄 LICENSE                          # License information
```

## 🔧 Requirements

- Python 3.7+
- Required packages (automatically installed):
  - `ete3` - Phylogenetic tree manipulation
  - `pandas` - Data analysis
  - `numpy` - Numerical computing
  - `scipy` - Statistical modeling
  - `matplotlib` - Basic plotting

## 📖 Documentation

- 📚 **[User Guide](docs/user_guide.md)** - Comprehensive usage instructions
- 🔬 **[Technical Documentation](docs/ChromEvol_Enhancement_README.md)** - Implementation details
- 📊 **[Scientific Report](docs/evolutionary_analysis_report.md)** - Analysis methodology and findings
- 🎓 **[Tutorial](examples/tutorial.md)** - Step-by-step examples

## 🎨 Example Output

The pipeline generates publication-ready visualizations with:
- Color-coded branches (blue=fusion, red=fission, green=stable)
- Statistical annotations (ML probabilities, support values)
- Multiple layout options (circular/rectangular)
- High-resolution output formats (SVG, PNG, PDF)

## 📊 Supported Analysis Types

### 1. Ancestral State Reconstruction
- Maximum likelihood estimation
- Parsimony-based inference
- Statistical confidence assessment

### 2. Chromosomal Rearrangement Analysis
- Fusion/fission event detection
- Event magnitude quantification
- Phylogenetic trend analysis

### 3. Model Selection & Optimization
- AIC/BIC model comparison
- Parameter optimization
- Rate heterogeneity modeling

### 4. Stochastic Mapping
- Event rate estimation
- Monte Carlo simulation
- Branch-specific statistics

## 🔬 Scientific Applications

This pipeline has been successfully applied to:
- Deuterostome chromosome evolution analysis
- Phylogenetic comparative genomics
- Karyotype evolution studies
- Genome organization research

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact & Support

- **Issues**: [GitHub Issues](https://github.com/your-username/chromevol-enhanced-analysis/issues)
- **Documentation**: [Project Wiki](https://github.com/your-username/chromevol-enhanced-analysis/wiki)
- **Email**: your.email@institution.edu

## 🏆 Citation

If you use this software in your research, please cite:

```bibtex
@software{chromevol_enhanced_2024,
  title={ChromEvol-Enhanced Chromosome Evolution Analysis},
  author={Your Name},
  year={2024},
  url={https://github.com/your-username/chromevol-enhanced-analysis},
  version={1.0.0}
}
```

## 🙏 Acknowledgments

- ChromEvol developers for methodological inspiration
- ETE3 team for phylogenetic tree tools
- Scientific Python community for foundational libraries

---

**⭐ If you find this project useful, please consider giving it a star!**
