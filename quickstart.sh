#!/bin/bash

# 🚀 ChromEvol-Enhanced Analysis Quick Start Script
# This script helps new users get started quickly with the analysis pipeline

set -e  # Exit on any error

echo "🧬 ChromEvol-Enhanced Chromosome Evolution Analysis"
echo "=================================================="
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is required but not found. Please install Python 3.7+ first."
    exit 1
fi

echo "✅ Python found: $(python3 --version)"

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "❌ pip3 is required but not found. Please install pip first."
    exit 1
fi

echo "✅ pip found"

# Install dependencies
echo ""
echo "📦 Installing required dependencies..."
pip3 install -r requirements.txt

# Check if data files exist
echo ""
echo "📁 Checking data files..."
if [ ! -f "data/pruned_phylogeny.nwk" ]; then
    echo "❌ Phylogenetic tree file not found: data/pruned_phylogeny.nwk"
    exit 1
fi

if [ ! -f "data/chromosome_counts.csv" ]; then
    echo "❌ Chromosome counts file not found: data/chromosome_counts.csv"
    exit 1
fi

echo "✅ Data files found"

# Create results directory if it doesn't exist
mkdir -p results

# Run basic analysis
echo ""
echo "🔬 Running basic chromosome evolution analysis..."
python3 src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --out_image results/quickstart_basic.svg \
    --out_ancestors results/quickstart_ancestors.csv \
    --out_rearrangements results/quickstart_events.csv

echo "✅ Basic analysis completed"

# Run ChromEvol analysis
echo ""
echo "🧬 Running ChromEvol-enhanced analysis..."
python3 src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --use_chromevol \
    --layout circular \
    --show_support \
    --out_image results/quickstart_chromevol.svg \
    --out_ancestors results/quickstart_chromevol_ancestors.csv \
    --out_rearrangements results/quickstart_chromevol_events.csv

echo "✅ ChromEvol analysis completed"

# Generate summary
echo ""
echo "📊 Generating analysis summary..."
python3 scripts/generate_summary.py

echo "✅ Summary generated"

# Display results
echo ""
echo "🎉 Analysis completed successfully!"
echo ""
echo "📋 Results generated:"
echo "   📊 Basic analysis: results/quickstart_basic.svg"
echo "   🧬 ChromEvol analysis: results/quickstart_chromevol.svg"
echo "   📈 Data summaries: results/analysis_summary_table.csv"
echo "   📝 Detailed reports: results/*_ancestors.csv, results/*_events.csv"
echo ""
echo "📖 Next steps:"
echo "   • View the generated SVG files to see your phylogenetic trees"
echo "   • Check the CSV reports for detailed numerical results"
echo "   • Read docs/user_guide.md for advanced usage options"
echo "   • Try examples/tutorial.md for a step-by-step walkthrough"
echo ""
echo "🔗 For help and documentation:"
echo "   • User guide: docs/user_guide.md"
echo "   • Tutorial: examples/tutorial.md"
echo "   • API reference: docs/api_reference.md"
echo "   • Issues: https://github.com/your-username/chromevol-enhanced-analysis/issues"
echo ""
echo "Happy analyzing! 🧬✨"
