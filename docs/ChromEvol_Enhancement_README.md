# ChromEvol-Enhanced Chromosome Evolution Analysis

## Overview

This enhanced version of the ancestral reconstruction script incorporates sophisticated methodologies inspired by ChromEvol, a state-of-the-art software for chromosome evolution analysis. The enhancement includes maximum likelihood methods, branch heterogeneity modeling, parametric optimization, and stochastic mapping.

## Key Features

### 1. ChromEvol-Inspired Maximum Likelihood Modeling

#### ChromEvolutionModel Class
- **Parametric evolution model** with rates for:
  - `gain`: Single chromosome gain events
  - `loss`: Single chromosome loss events  
  - `dupl`: Whole genome duplication events
  - `demidupl`: Polyploidization followed by diploidization
  - `linear_rate`: Chromosome number-dependent rate modification

#### Transition Rate Matrix
- **Q-matrix construction** following ChromEvol principles
- **Matrix exponentiation** for transition probabilities: P(t) = exp(Q*t)
- **Branch-specific parameters** for heterogeneous models

### 2. Maximum Likelihood Optimization

#### ChromEvolOptimizer Class
- **Parameter optimization** using scipy.optimize
- **Multiple optimization methods** (L-BFGS-B, etc.)
- **Bounded optimization** ensuring positive rates
- **Convergence diagnostics** and success reporting

#### Model Selection
- **AIC/BIC criteria** for model comparison
- **Multiple model types**: simple vs. full models
- **Automatic best model selection** based on information criteria

### 3. Likelihood-Based Ancestral State Reconstruction

#### ChromEvolLikelihoodCalculator Class
- **Felsenstein's pruning algorithm** for likelihood calculation
- **Marginal ancestral state reconstruction** 
- **Maximum likelihood state assignment** with probabilities
- **Root frequency constraints** (can specify prior beliefs)

### 4. Stochastic Mapping

#### StochasticMapper Class
- **Event simulation** along branches using rate matrices
- **Monte Carlo sampling** (default: 100 replicates)
- **Event type classification**: gains, losses, duplications, other
- **Confidence intervals** for expected event counts
- **Branch-specific event summaries**

## Command Line Interface

### Basic Usage

```bash
# Traditional parsimony-based analysis (default)
python ancestral_reconstruction.py

# ChromEvol-inspired maximum likelihood analysis
python ancestral_reconstruction.py --use_chromevol

# Full optimization with model selection
python ancestral_reconstruction.py --use_chromevol --optimize_params --model_selection
```

### New Parameters

#### ChromEvol Features
- `--use_chromevol`: Enable maximum likelihood methods instead of parsimony
- `--optimize_params`: Optimize model parameters using maximum likelihood
- `--model_selection`: Perform AIC/BIC-based model selection

#### Existing Parameters (Enhanced)
- `--root_count`: Constrains root chromosome number (used as prior in ML)
- `--out_ancestors`: Now includes ML probabilities and expected events
- `--out_image`: Enhanced visualization with event information

## Output Enhancements

### Ancestral States Report
Enhanced CSV output now includes:
- `node_name`: Internal node identifier
- `inferred_chromosome_count`: Most likely ancestral state
- `ml_probability`: Maximum likelihood probability of the state
- `expected_gains`: Expected number of gain events (from stochastic mapping)
- `expected_losses`: Expected number of loss events
- `expected_duplications`: Expected number of duplication events

### Enhanced Visualization
- **Maximum likelihood probabilities** displayed on internal nodes
- **Expected event counts** shown as branch annotations
- **Orange highlighting** for branches with high duplication activity
- **Enhanced legend** explaining all visual elements

### Console Output
- **Parameter optimization progress** with initial and final values
- **Log-likelihood values** for model comparison
- **AIC/BIC model selection** results
- **Stochastic mapping summaries** with expected event counts per branch

## Scientific Background

### ChromEvol Inspiration
This implementation draws from the ChromEvol methodology (Mayrose et al. 2010, Glick & Mayrose 2014):

1. **Maximum Likelihood Framework**: Uses probabilistic models instead of parsimony
2. **Parametric Models**: Explicit rates for different types of chromosome changes
3. **Branch Heterogeneity**: Allows different parameters on different branches
4. **Stochastic Mapping**: Simulates evolutionary histories to estimate event expectations

### Mathematical Framework

#### Rate Matrix Construction
The instantaneous rate matrix Q follows ChromEvol's approach:
- Q[i,i+1] = gain_rate * (1 + linear_rate * i)
- Q[i,i-1] = loss_rate * (1 + linear_rate * i)  
- Q[i,2*i] = duplication_rate
- Q[i,j] += demiduplication_rate (for nearby states)

#### Likelihood Calculation
Total likelihood: L = Σ P(root=i) * L(data|root=i)
Where L(data|root=i) is calculated using Felsenstein's algorithm.

#### Optimization
Parameters are optimized to maximize log-likelihood using gradient-based methods with bounds.

## Performance Considerations

### Computational Complexity
- **Matrix exponentiation**: O(n³) where n is max chromosome number
- **Stochastic mapping**: O(replicates × branches × events)
- **Optimization**: Depends on convergence and parameter space

### Memory Usage
- Transition matrices: O(n²) per branch
- Likelihood arrays: O(nodes × n)
- Stochastic mapping: O(replicates × branches × events)

## Validation and Testing

### Test Cases
The enhanced script has been tested with:
- Small trees (5-10 species)
- Medium trees (20-50 species) 
- Large trees (100+ species)
- Various chromosome number ranges (5-50)

### Comparison with ChromEvol
Key differences from original ChromEvol:
- **Simplified transition model**: Focuses on core events
- **Python implementation**: More accessible than C++
- **Integrated visualization**: Built-in tree plotting
- **Streamlined interface**: Single script execution

## Future Enhancements

### Potential Extensions
1. **Base number optimization**: Automatic inference of base chromosome number
2. **Rate heterogeneity**: Gamma-distributed rates across branches
3. **Model complexity**: Additional transition types (translocations, etc.)
4. **Bayesian inference**: MCMC sampling for parameter uncertainty
5. **Comparative analysis**: Multi-gene family chromosome evolution

### Integration Opportunities
- **R integration**: Export for detailed statistical analysis
- **Visualization enhancement**: Interactive plots with plotly
- **Database connectivity**: Direct access to chromosome databases
- **Pipeline integration**: Snakemake/Nextflow workflow compatibility

## References

1. Mayrose, I., Barker, M. S., & Otto, S. P. (2010). Probabilistic models of chromosome number evolution and the inference of polyploidy. *Systematic Biology*, 59(2), 132-144.

2. Glick, L., & Mayrose, I. (2014). ChromEvol: assessing the pattern of chromosome number evolution and the inference of polyploidy along a phylogeny. *Molecular Biology and Evolution*, 31(7), 1914-1922.

3. Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17(6), 368-376.

## Citation

If you use this enhanced script in your research, please cite:
- The original ChromEvol papers (Mayrose et al. 2010, Glick & Mayrose 2014)
- Your usage of this enhanced implementation

## Support

For questions, bug reports, or feature requests related to the ChromEvol enhancements, please refer to the script's inline documentation and comments.
