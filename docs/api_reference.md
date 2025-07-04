# 🔧 API Reference

## 📚 Module Overview

The ChromEvol-Enhanced Analysis toolkit consists of several core modules and classes for chromosome evolution analysis.

## 🧬 Core Classes

### ChromEvolutionModel

Maximum likelihood chromosome evolution model inspired by ChromEvol methodology.

```python
class ChromEvolutionModel:
    def __init__(self, max_chromosome_number=50, model_type='full')
```

**Parameters:**
- `max_chromosome_number` (int): Maximum chromosome number to consider in the model
- `model_type` (str): Model complexity ('simple', 'full', 'gain_loss_only')

**Methods:**

#### set_parameters(**kwargs)
Set model evolutionary rate parameters.

```python
model.set_parameters(gain=0.1, loss=0.2, dupl=0.05)
```

**Parameters:**
- `gain` (float): Rate of single chromosome gain events
- `loss` (float): Rate of single chromosome loss events  
- `dupl` (float): Rate of whole genome duplication events
- `demidupl` (float): Rate of polyploidization followed by diploidization
- `base_number` (int): Base chromosome number for the model
- `linear_rate` (float): Linear chromosome-number dependent rate

#### build_rate_matrix(chromosome_range)
Construct the Q-matrix for transition rates between chromosome states.

```python
Q_matrix = model.build_rate_matrix(range(1, 51))
```

**Parameters:**
- `chromosome_range` (range): Range of chromosome numbers to include

**Returns:**
- `numpy.ndarray`: Transition rate matrix (Q-matrix)

#### calculate_transition_probabilities(branch_length)
Calculate transition probabilities using matrix exponentiation.

```python
P_matrix = model.calculate_transition_probabilities(0.1)
```

**Parameters:**
- `branch_length` (float): Phylogenetic branch length

**Returns:**
- `numpy.ndarray`: Transition probability matrix P(t) = exp(Q×t)

### ChromEvolLikelihoodCalculator

Implements Felsenstein's pruning algorithm for likelihood calculation.

```python
class ChromEvolLikelihoodCalculator:
    def __init__(self, model, tree, chromosome_data)
```

**Parameters:**
- `model` (ChromEvolutionModel): Evolution model instance
- `tree` (ete3.Tree): Phylogenetic tree
- `chromosome_data` (dict): Species chromosome count data

**Methods:**

#### calculate_likelihood(parameters=None)
Calculate the likelihood of the data given model parameters.

```python
log_likelihood = calculator.calculate_likelihood()
```

**Parameters:**
- `parameters` (dict, optional): Model parameters to use

**Returns:**
- `float`: Log-likelihood value

#### reconstruct_ancestral_states()
Perform marginal ancestral state reconstruction.

```python
ancestral_states = calculator.reconstruct_ancestral_states()
```

**Returns:**
- `dict`: Dictionary mapping node names to inferred chromosome counts

#### get_state_probabilities(node_name)
Get the probability distribution over chromosome states for a node.

```python
probabilities = calculator.get_state_probabilities("internal_node_1")
```

**Parameters:**
- `node_name` (str): Name of the phylogenetic node

**Returns:**
- `numpy.ndarray`: Probability distribution over chromosome states

### ChromEvolOptimizer

Optimization engine for maximum likelihood parameter estimation.

```python
class ChromEvolOptimizer:
    def __init__(self, likelihood_calculator)
```

**Parameters:**
- `likelihood_calculator` (ChromEvolLikelihoodCalculator): Likelihood calculator instance

**Methods:**

#### optimize_parameters(initial_params=None, bounds=None, method='L-BFGS-B')
Optimize model parameters using maximum likelihood.

```python
result = optimizer.optimize_parameters(
    initial_params={'gain': 0.1, 'loss': 0.1},
    bounds={'gain': (0, 1), 'loss': (0, 1)}
)
```

**Parameters:**
- `initial_params` (dict): Starting parameter values
- `bounds` (dict): Parameter bounds for optimization
- `method` (str): Optimization algorithm to use

**Returns:**
- `OptimizationResult`: Optimization result object with optimized parameters

#### calculate_model_selection_criteria()
Calculate AIC and BIC for model selection.

```python
aic, bic = optimizer.calculate_model_selection_criteria()
```

**Returns:**
- `tuple`: (AIC value, BIC value)

### StochasticMapper

Performs stochastic mapping to estimate expected numbers of evolutionary events.

```python
class StochasticMapper:
    def __init__(self, model, tree, ancestral_states)
```

**Parameters:**
- `model` (ChromEvolutionModel): Evolution model instance
- `tree` (ete3.Tree): Phylogenetic tree
- `ancestral_states` (dict): Reconstructed ancestral states

**Methods:**

#### perform_stochastic_mapping(num_replicates=100)
Simulate evolutionary events along branches.

```python
event_counts = mapper.perform_stochastic_mapping(num_replicates=1000)
```

**Parameters:**
- `num_replicates` (int): Number of Monte Carlo replicates

**Returns:**
- `dict`: Expected event counts per branch

#### simulate_branch_events(start_state, end_state, branch_length)
Simulate events along a single branch.

```python
events = mapper.simulate_branch_events(start_state=20, end_state=22, branch_length=0.1)
```

**Parameters:**
- `start_state` (int): Starting chromosome number
- `end_state` (int): Ending chromosome number  
- `branch_length` (float): Phylogenetic branch length

**Returns:**
- `dict`: Dictionary of event counts by type

## 🔧 Utility Functions

### Data Loading and Validation

#### load_phylogenetic_tree(file_path)
Load and validate phylogenetic tree from Newick format.

```python
from src.ancestral_reconstruction import load_phylogenetic_tree
tree = load_phylogenetic_tree("data/tree.nwk")
```

**Parameters:**
- `file_path` (str): Path to Newick format tree file

**Returns:**
- `ete3.Tree`: Loaded phylogenetic tree

#### load_chromosome_data(file_path)
Load chromosome count data from CSV file.

```python
from src.ancestral_reconstruction import load_chromosome_data
data = load_chromosome_data("data/counts.csv")
```

**Parameters:**
- `file_path` (str): Path to CSV file with species and chromosome counts

**Returns:**
- `dict`: Dictionary mapping species names to chromosome counts

#### validate_data_consistency(tree, chromosome_data)
Check consistency between tree and chromosome data.

```python
from src.ancestral_reconstruction import validate_data_consistency
is_valid = validate_data_consistency(tree, chromosome_data)
```

**Parameters:**
- `tree` (ete3.Tree): Phylogenetic tree
- `chromosome_data` (dict): Chromosome count data

**Returns:**
- `bool`: True if data is consistent

### Analysis Functions

#### perform_parsimony_reconstruction(tree, chromosome_data)
Perform maximum parsimony ancestral state reconstruction.

```python
from src.ancestral_reconstruction import perform_parsimony_reconstruction
ancestral_states = perform_parsimony_reconstruction(tree, chromosome_data)
```

**Parameters:**
- `tree` (ete3.Tree): Phylogenetic tree
- `chromosome_data` (dict): Chromosome count data

**Returns:**
- `dict`: Inferred ancestral chromosome states

#### perform_chromevol_analysis(tree, chromosome_data, **kwargs)
Perform complete ChromEvol-style analysis.

```python
from src.ancestral_reconstruction import perform_chromevol_analysis
results = perform_chromevol_analysis(
    tree, 
    chromosome_data,
    optimize_params=True,
    model_selection=True
)
```

**Parameters:**
- `tree` (ete3.Tree): Phylogenetic tree
- `chromosome_data` (dict): Chromosome count data
- `optimize_params` (bool): Whether to optimize parameters
- `model_selection` (bool): Whether to perform model selection

**Returns:**
- `dict`: Complete analysis results

#### analyze_rearrangement_events(tree, ancestral_states, chromosome_data)
Quantify chromosomal rearrangement events.

```python
from src.ancestral_reconstruction import analyze_rearrangement_events
events = analyze_rearrangement_events(tree, ancestral_states, chromosome_data)
```

**Parameters:**
- `tree` (ete3.Tree): Phylogenetic tree
- `ancestral_states` (dict): Reconstructed ancestral states
- `chromosome_data` (dict): Terminal chromosome counts

**Returns:**
- `list`: List of rearrangement events with details

### Visualization Functions

#### plot_annotated_tree(tree, ancestral_states, **kwargs)
Generate annotated phylogenetic tree visualization.

```python
from src.ancestral_reconstruction import plot_annotated_tree
plot_annotated_tree(
    tree, 
    ancestral_states,
    layout='circular',
    show_support=True,
    show_branch_length=True,
    output_file='tree.svg'
)
```

**Parameters:**
- `tree` (ete3.Tree): Phylogenetic tree
- `ancestral_states` (dict): Ancestral chromosome states
- `layout` (str): Tree layout ('circular' or 'rectangular')
- `show_support` (bool): Display branch support values
- `show_branch_length` (bool): Display branch lengths
- `output_file` (str): Output file path

**Returns:**
- `None`: Saves visualization to file

## 📊 Data Structures

### Optimization Result

Result object returned by parameter optimization.

```python
class OptimizationResult:
    success: bool              # Whether optimization succeeded
    parameters: dict           # Optimized parameter values  
    log_likelihood: float      # Final log-likelihood
    iterations: int            # Number of optimization iterations
    message: str              # Optimization status message
```

### Analysis Configuration

Configuration object for analysis parameters.

```python
class AnalysisConfig:
    use_chromevol: bool        # Use ChromEvol methods
    optimize_params: bool      # Optimize parameters
    model_selection: bool      # Perform model selection
    root_constraint: int       # Fixed root chromosome number
    max_iterations: int        # Maximum optimization iterations
    convergence_tol: float     # Convergence tolerance
```

## 🎨 Visualization Parameters

### Tree Styling Options

```python
tree_style_params = {
    'layout': 'circular',           # 'circular' or 'rectangular'
    'show_support': True,           # Display support values
    'show_branch_length': True,     # Display branch lengths
    'node_size': 10,               # Node circle size
    'branch_width': 2,             # Branch line width
    'font_size': 12,               # Text font size
    'dpi': 300,                    # Image resolution
    'width': 800,                  # Image width (pixels)
    'height': 600                  # Image height (pixels)
}
```

### Color Schemes

```python
color_scheme = {
    'fusion': '#0066CC',        # Blue for fusion events
    'fission': '#CC0000',       # Red for fission events  
    'stable': '#00AA00',        # Green for stable branches
    'duplication': '#FF8800',   # Orange for high duplication
    'uncertain': '#CCCCCC'      # Gray for uncertain data
}
```

## ⚠️ Error Handling

### Common Exceptions

```python
class ChromEvolError(Exception):
    """Base exception for ChromEvol analysis errors"""
    pass

class OptimizationError(ChromEvolError):
    """Raised when parameter optimization fails"""
    pass

class DataValidationError(ChromEvolError):
    """Raised when input data validation fails"""
    pass

class TreeFormatError(ChromEvolError):
    """Raised when phylogenetic tree format is invalid"""
    pass
```

### Error Handling Example

```python
try:
    result = perform_chromevol_analysis(tree, data, optimize_params=True)
except OptimizationError as e:
    print(f"Optimization failed: {e}")
    # Fall back to non-optimized analysis
    result = perform_chromevol_analysis(tree, data, optimize_params=False)
except DataValidationError as e:
    print(f"Data validation error: {e}")
    # Handle data inconsistency
```

## 📈 Performance Notes

### Memory Usage
- Memory usage scales with O(n × m) where n = number of species, m = chromosome range
- Large trees (>100 species) may require 8GB+ RAM for optimization
- Use `max_chromosome_number` parameter to limit memory usage

### Computational Complexity
- Parsimony reconstruction: O(n × m)
- Maximum likelihood: O(n × m² × iterations)
- Stochastic mapping: O(replicates × n × m)

### Optimization Tips
- Start with parsimony for initial estimates
- Use smaller chromosome ranges for large datasets
- Enable parallel processing for stochastic mapping

---

For more detailed examples and usage patterns, see the [User Guide](user_guide.md) and [Tutorial](../examples/tutorial.md).
