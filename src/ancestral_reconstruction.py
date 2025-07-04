import pandas as pd
from ete3 import Tree
from ete3.treeview import TreeStyle, NodeStyle, TextFace
import argparse
import os
import numpy as np
import scipy.optimize as opt
import scipy.linalg
from scipy.special import logsumexp
from scipy.stats import poisson
import itertools
import warnings
warnings.filterwarnings("ignore")

# ChromEvol-inspired classes and functions
class ChromEvolutionModel:
    """
    ChromEvol-inspired chromosome evolution model with maximum likelihood estimation,
    branch heterogeneity, and parametric modeling of gain, loss, duplication events.
    """
    
    def __init__(self, max_chromosome_number=50, model_type='full'):
        self.max_chr = max_chromosome_number
        self.model_type = model_type  # 'full', 'simple', 'gain_loss_only'
        
        # Model parameters (inspired by ChromEvol)
        self.params = {
            'gain': 0.1,      # Rate of chromosome gain
            'loss': 0.1,      # Rate of chromosome loss  
            'dupl': 0.05,     # Rate of whole genome duplication
            'demidupl': 0.02, # Rate of polyploidization followed by diploidization
            'base_number': 7, # Base chromosome number
            'linear_rate': 0.01  # Linear chromosome-number dependency
        }
        
        # Branch-specific parameters for heterogeneity
        self.branch_params = {}
        
        # Transition rate matrix
        self.Q_matrix = None
        self.transition_probs = {}
        
    def set_parameters(self, **kwargs):
        """Update model parameters"""
        for key, value in kwargs.items():
            if key in self.params:
                self.params[key] = value
    
    def set_branch_parameters(self, branch_id, **kwargs):
        """Set branch-specific parameters for heterogeneous models"""
        if branch_id not in self.branch_params:
            self.branch_params[branch_id] = self.params.copy()
        for key, value in kwargs.items():
            if key in self.branch_params[branch_id]:
                self.branch_params[branch_id][key] = value
    
    def build_transition_rate_matrix(self, params=None):
        """
        Build the instantaneous rate matrix Q for chromosome number transitions
        Based on ChromEvol's approach with gain, loss, duplication events
        """
        if params is None:
            params = self.params
            
        Q = np.zeros((self.max_chr + 1, self.max_chr + 1))
        
        for i in range(self.max_chr + 1):
            # Chromosome gain (i -> i+1)
            if i < self.max_chr:
                rate = params['gain']
                if 'linear_rate' in params:
                    rate *= (1 + params['linear_rate'] * i)  # Linear dependency
                Q[i, i + 1] = rate
            
            # Chromosome loss (i -> i-1) 
            if i > 1:  # Cannot go below 1 chromosome
                rate = params['loss']
                if 'linear_rate' in params:
                    rate *= (1 + params['linear_rate'] * i)  # Linear dependency
                Q[i, i - 1] = rate
            
            # Whole genome duplication (i -> 2*i)
            if 2 * i <= self.max_chr and i > 0:
                Q[i, 2 * i] = params['dupl']
            
            # Demi-duplication (polyploidization followed by loss)
            if 'demidupl' in params:
                # Simplified: add rate to nearby states
                for j in [i-2, i-1, i+1, i+2]:
                    if 0 < j <= self.max_chr and j != i:
                        Q[i, j] += params['demidupl'] / 4
        
        # Set diagonal elements to make rows sum to zero
        for i in range(self.max_chr + 1):
            Q[i, i] = -np.sum(Q[i, :])
            
        self.Q_matrix = Q
        return Q
    
    def get_transition_probabilities(self, branch_length, params=None):
        """
        Calculate transition probabilities P(t) = exp(Q*t) for a given branch length
        """
        if params is None:
            Q = self.build_transition_rate_matrix()
        else:
            Q = self.build_transition_rate_matrix(params)
            
        try:
            # Matrix exponentiation: P(t) = exp(Q*t)
            P = scipy.linalg.expm(Q * branch_length)
            # Ensure probabilities are non-negative and normalized
            P = np.maximum(P, 1e-10)
            P = P / P.sum(axis=1, keepdims=True)
            return P
        except:
            # Fallback to simple exponential decay if matrix exponentiation fails
            P = np.eye(self.max_chr + 1)
            for i in range(self.max_chr + 1):
                for j in range(self.max_chr + 1):
                    if i != j:
                        P[i, j] = np.exp(-branch_length) * 0.01
                P[i, i] = 1 - np.sum(P[i, :i]) - np.sum(P[i, i+1:])
            return P

class ChromEvolLikelihoodCalculator:
    """
    Maximum likelihood calculator for chromosome evolution following ChromEvol principles
    """
    
    def __init__(self, tree, chromosome_data, model):
        self.tree = tree
        self.data = chromosome_data
        self.model = model
        
        # Likelihood matrices for each node
        self.likelihoods = {}
        
        # Root frequencies
        self.root_frequencies = None
        
    def set_root_frequencies(self, frequencies):
        """Set root state frequencies"""
        self.root_frequencies = frequencies
    
    def calculate_likelihood_at_node(self, node):
        """
        Calculate likelihood at a given node using Felsenstein's pruning algorithm
        """
        if node.is_leaf():
            # For leaves, likelihood is 1 for observed state, 0 for others
            likelihood = np.zeros(self.model.max_chr + 1)
            if hasattr(node, 'count') and node.count is not None:
                if 0 <= node.count <= self.model.max_chr:
                    likelihood[node.count] = 1.0
                else:
                    # If observed count is outside our range, distribute probability
                    likelihood.fill(1.0 / len(likelihood))
            else:
                # Missing data - equal probability for all states
                likelihood.fill(1.0 / len(likelihood))
            
            self.likelihoods[id(node)] = likelihood
            return likelihood
        
        # For internal nodes, multiply child likelihoods with transition probabilities
        likelihood = np.ones(self.model.max_chr + 1)
        
        for child in node.children:
            child_likelihood = self.calculate_likelihood_at_node(child)
            
            # Get transition probabilities for this branch
            branch_length = child.dist if hasattr(child, 'dist') else 0.1
            
            # Use branch-specific parameters if available
            branch_id = id(child)
            if branch_id in self.model.branch_params:
                P = self.model.get_transition_probabilities(branch_length, 
                                                          self.model.branch_params[branch_id])
            else:
                P = self.model.get_transition_probabilities(branch_length)
            
            # Calculate conditional likelihood
            conditional_likelihood = np.dot(P, child_likelihood)
            likelihood *= conditional_likelihood
        
        self.likelihoods[id(node)] = likelihood
        return likelihood
    
    def calculate_total_likelihood(self):
        """
        Calculate total likelihood of the tree
        """
        root = self.tree.get_tree_root()
        root_likelihood = self.calculate_likelihood_at_node(root)
        
        if self.root_frequencies is None:
            # Uniform root frequencies if not specified
            self.root_frequencies = np.ones(self.model.max_chr + 1) / (self.model.max_chr + 1)
        
        # Total likelihood is sum over all possible root states
        total_likelihood = np.sum(root_likelihood * self.root_frequencies)
        return np.log(total_likelihood) if total_likelihood > 0 else -np.inf
    
    def ancestral_state_reconstruction(self):
        """
        Perform marginal ancestral state reconstruction (maximum likelihood states)
        """
        # First calculate likelihoods
        root = self.tree.get_tree_root()
        self.calculate_likelihood_at_node(root)
        
        # Reconstruct ancestral states using marginal reconstruction
        for node in self.tree.traverse("postorder"):
            if not node.is_leaf():
                likelihood = self.likelihoods[id(node)]
                # Find most likely state
                most_likely_state = np.argmax(likelihood)
                node.add_features(ml_count=most_likely_state)
                node.add_features(ml_probability=likelihood[most_likely_state])

class ChromEvolOptimizer:
    """
    Parameter optimization for chromosome evolution models using maximum likelihood
    """
    
    def __init__(self, tree, chromosome_data, model):
        self.tree = tree
        self.data = chromosome_data
        self.model = model
        self.calculator = ChromEvolLikelihoodCalculator(tree, chromosome_data, model)
    
    def objective_function(self, params_vector):
        """
        Objective function for optimization (negative log-likelihood)
        """
        try:
            # Update model parameters
            param_names = ['gain', 'loss', 'dupl', 'demidupl', 'linear_rate']
            for i, param_name in enumerate(param_names):
                if i < len(params_vector):
                    self.model.set_parameters(**{param_name: max(params_vector[i], 1e-6)})
            
            # Calculate likelihood
            log_likelihood = self.calculator.calculate_total_likelihood()
            return -log_likelihood  # Return negative for minimization
        except:
            return 1e6  # Return large value if calculation fails
    
    def optimize_parameters(self, method='L-BFGS-B', maxiter=1000):
        """
        Optimize model parameters using maximum likelihood
        """
        # Initial parameter values
        initial_params = [
            self.model.params['gain'],
            self.model.params['loss'], 
            self.model.params['dupl'],
            self.model.params['demidupl'],
            self.model.params['linear_rate']
        ]
        
        # Parameter bounds (all rates must be positive)
        bounds = [(1e-6, 10)] * len(initial_params)
        
        print("Optimizing model parameters using maximum likelihood...")
        print(f"Initial parameters: gain={initial_params[0]:.4f}, loss={initial_params[1]:.4f}, "
              f"dupl={initial_params[2]:.4f}, demidupl={initial_params[3]:.4f}, linear={initial_params[4]:.4f}")
        
        # Optimize
        result = opt.minimize(
            self.objective_function,
            initial_params,
            method=method,
            bounds=bounds,
            options={'maxiter': maxiter}
        )
        
        if result.success:
            # Update model with optimized parameters
            param_names = ['gain', 'loss', 'dupl', 'demidupl', 'linear_rate']
            optimized_params = {}
            for i, param_name in enumerate(param_names):
                if i < len(result.x):
                    optimized_params[param_name] = result.x[i]
            
            self.model.set_parameters(**optimized_params)
            
            print(f"Optimization successful!")
            print(f"Optimized parameters: gain={result.x[0]:.4f}, loss={result.x[1]:.4f}, "
                  f"dupl={result.x[2]:.4f}, demidupl={result.x[3]:.4f}, linear={result.x[4]:.4f}")
            print(f"Final log-likelihood: {-result.fun:.2f}")
            
            return result
        else:
            print(f"Optimization failed: {result.message}")
            return None

class StochasticMapper:
    """
    Stochastic mapping implementation for chromosome evolution events
    """
    
    def __init__(self, tree, model, num_mappings=100):
        self.tree = tree
        self.model = model
        self.num_mappings = num_mappings
        self.event_counts = {}
    
    def simulate_events_on_branch(self, start_state, end_state, branch_length, params):
        """
        Simulate chromosome evolution events along a branch using stochastic mapping
        """
        events = []
        current_state = start_state
        current_time = 0.0
        
        # Build rate matrix for this branch
        Q = self.model.build_transition_rate_matrix(params)
        
        while current_time < branch_length and current_state != end_state:
            # Calculate total rate out of current state
            total_rate = -Q[current_state, current_state]
            if total_rate <= 0:
                break
                
            # Sample waiting time
            waiting_time = np.random.exponential(1.0 / total_rate)
            if current_time + waiting_time > branch_length:
                break
                
            current_time += waiting_time
            
            # Sample next state
            rates = Q[current_state, :].copy()
            rates[current_state] = 0  # Remove diagonal
            rates = np.maximum(rates, 0)  # Ensure non-negative
            
            if np.sum(rates) > 0:
                probs = rates / np.sum(rates)
                next_state = np.random.choice(len(probs), p=probs)
                
                # Classify event type
                if next_state == current_state + 1:
                    event_type = 'gain'
                elif next_state == current_state - 1:
                    event_type = 'loss'
                elif next_state == 2 * current_state:
                    event_type = 'duplication'
                else:
                    event_type = 'other'
                
                events.append({
                    'time': current_time,
                    'from_state': current_state,
                    'to_state': next_state,
                    'type': event_type
                })
                
                current_state = next_state
            else:
                break
        
        return events
    
    def perform_stochastic_mapping(self):
        """
        Perform stochastic mapping across the entire tree
        """
        all_mappings = []
        
        print(f"Performing stochastic mapping with {self.num_mappings} replicates...")
        
        for mapping_i in range(self.num_mappings):
            mapping_events = {}
            
            for node in self.tree.traverse():
                if not node.is_root() and hasattr(node, 'ml_count') and hasattr(node.up, 'ml_count'):
                    branch_length = node.dist if hasattr(node, 'dist') else 0.1
                    branch_id = id(node)
                    
                    # Use branch-specific parameters if available
                    params = self.model.branch_params.get(branch_id, self.model.params)
                    
                    # Simulate events on this branch
                    events = self.simulate_events_on_branch(
                        node.up.ml_count, node.ml_count, branch_length, params
                    )
                    
                    mapping_events[branch_id] = events
            
            all_mappings.append(mapping_events)
        
        # Summarize event counts
        self.summarize_event_counts(all_mappings)
        return all_mappings
    
    def summarize_event_counts(self, all_mappings):
        """
        Summarize event counts across all stochastic mappings
        """
        event_summary = {}
        
        for node in self.tree.traverse():
            if not node.is_root():
                branch_id = id(node)
                event_types = ['gain', 'loss', 'duplication', 'other']
                event_counts = {event_type: [] for event_type in event_types}
                
                for mapping in all_mappings:
                    if branch_id in mapping:
                        branch_events = mapping[branch_id]
                        for event_type in event_types:
                            count = sum(1 for event in branch_events if event['type'] == event_type)
                            event_counts[event_type].append(count)
                    else:
                        for event_type in event_types:
                            event_counts[event_type].append(0)
                
                # Calculate means and confidence intervals
                summary = {}
                for event_type in event_types:
                    counts = event_counts[event_type]
                    summary[event_type] = {
                        'mean': np.mean(counts),
                        'std': np.std(counts),
                        'ci_lower': np.percentile(counts, 2.5),
                        'ci_upper': np.percentile(counts, 97.5)
                    }
                
                event_summary[branch_id] = summary
                
                # Add to node for visualization
                node.add_features(expected_gains=summary['gain']['mean'])
                node.add_features(expected_losses=summary['loss']['mean'])
                node.add_features(expected_duplications=summary['duplication']['mean'])
        
        self.event_counts = event_summary
        return event_summary

# Enhanced analysis functions

def perform_chromevol_analysis(tree, counts_dict, use_ml=True, optimize_params=True, stochastic_mapping=True):
    """
    Perform ChromEvol-inspired maximum likelihood analysis of chromosome evolution
    """
    print("\n--- ChromEvol-Inspired Maximum Likelihood Analysis ---")
    
    # Initialize model
    model = ChromEvolutionModel(max_chromosome_number=50, model_type='full')
    
    # Attach data to tree
    for leaf in tree:
        if leaf.name in counts_dict:
            leaf.add_features(count=counts_dict[leaf.name])
        else:
            leaf.add_features(count=None)
    
    if use_ml:
        # Maximum likelihood optimization
        optimizer = ChromEvolOptimizer(tree, counts_dict, model)
        
        if optimize_params:
            # Optimize parameters
            result = optimizer.optimize_parameters()
        
        # Perform ML ancestral state reconstruction
        calculator = ChromEvolLikelihoodCalculator(tree, counts_dict, model)
        
        # Set root frequencies (can be optimized or set based on prior knowledge)
        root_freq = np.zeros(model.max_chr + 1)
        root_freq[24] = 1.0  # Strong prior on root state = 24
        calculator.set_root_frequencies(root_freq)
        
        # Calculate final likelihood
        final_likelihood = calculator.calculate_total_likelihood()
        print(f"Final log-likelihood: {final_likelihood:.2f}")
        
        # Perform ancestral reconstruction
        calculator.ancestral_state_reconstruction()
        
        # Copy ML results to standard count attribute for visualization
        for node in tree.traverse():
            if hasattr(node, 'ml_count'):
                node.add_features(count=node.ml_count)
                node.add_features(ml_probability=getattr(node, 'ml_probability', 0.0))
        
        if stochastic_mapping:
            # Perform stochastic mapping
            mapper = StochasticMapper(tree, model, num_mappings=100)
            mappings = mapper.perform_stochastic_mapping()
            
            # Generate stochastic mapping report
            print("\n--- Expected Number of Events per Branch (Stochastic Mapping) ---")
            for node in tree.traverse():
                if not node.is_root() and hasattr(node, 'expected_gains'):
                    node_name = node.name if node.name else f"Branch_to_{node.get_leaves()[0].name}"
                    print(f"{node_name}: Gains={node.expected_gains:.2f}, "
                          f"Losses={node.expected_losses:.2f}, "
                          f"Duplications={node.expected_duplications:.2f}")
    
    return model

def calculate_model_selection_criteria(tree, counts_dict, models=['simple', 'full']):
    """
    Calculate AIC/BIC for model selection
    """
    print("\n--- Model Selection ---")
    
    results = {}
    
    for model_type in models:
        print(f"Evaluating {model_type} model...")
        
        model = ChromEvolutionModel(max_chromosome_number=50, model_type=model_type)
        optimizer = ChromEvolOptimizer(tree, counts_dict, model)
        
        # Optimize parameters
        opt_result = optimizer.optimize_parameters(maxiter=500)
        
        if opt_result and opt_result.success:
            log_likelihood = -opt_result.fun
            n_params = len(opt_result.x)
            n_species = len([n for n in tree if n.is_leaf()])
            
            # Calculate AIC and BIC
            AIC = 2 * n_params - 2 * log_likelihood
            BIC = np.log(n_species) * n_params - 2 * log_likelihood
            
            results[model_type] = {
                'log_likelihood': log_likelihood,
                'AIC': AIC,
                'BIC': BIC,
                'n_params': n_params,
                'parameters': dict(zip(['gain', 'loss', 'dupl', 'demidupl', 'linear_rate'], opt_result.x))
            }
            
            print(f"  Log-likelihood: {log_likelihood:.2f}")
            print(f"  AIC: {AIC:.2f}")
            print(f"  BIC: {BIC:.2f}")
        else:
            print(f"  Optimization failed for {model_type} model")
    
    # Select best model
    if results:
        best_aic = min(results.keys(), key=lambda x: results[x]['AIC'])
        best_bic = min(results.keys(), key=lambda x: results[x]['BIC'])
        
        print(f"\nBest model by AIC: {best_aic}")
        print(f"Best model by BIC: {best_bic}")
    
    return results

def plot_annotated_tree(tree, output_file, layout, show_branch_length, show_support, dpi):
    """
    Generates and saves a high-quality, aesthetically pleasing image of the 
    phylogenetic tree with annotations.
    """
    print("\n--- Generating Annotated Phylogenetic Tree Image ---")
    try:
        # --- Create Tree Style ---
        ts = TreeStyle()
        ts.mode = "c" if layout == "circular" else "r"  # Circular or rectangular layout
        ts.show_leaf_name = True
        ts.show_branch_length = show_branch_length  # User-configurable
        ts.show_branch_support = show_support  # User-configurable branch support display
        ts.scale = 120

        # --- Add Title ---
        ts.title.add_face(TextFace("Ancestral Chromosome Evolution in Echinoderms", fsize=20, fgcolor="black", bold=True), column=0)

        # --- Add Enhanced Legend for ChromEvol Features ---
        ts.legend.add_face(TextFace("■", fsize=15, fgcolor="#377eb8"), column=0)
        ts.legend.add_face(TextFace(" Fusion Event (n decreases)", fsize=12), column=1)
        ts.legend.add_face(TextFace("■", fsize=15, fgcolor="#e41a1c"), column=0)
        ts.legend.add_face(TextFace(" Fission Event (n increases)", fsize=12), column=1)
        ts.legend.add_face(TextFace("■", fsize=15, fgcolor="#4daf4a"), column=0)
        ts.legend.add_face(TextFace(" Stable Lineage", fsize=12), column=1)
        ts.legend.add_face(TextFace("■", fsize=15, fgcolor="#ff7f00"), column=0)
        ts.legend.add_face(TextFace(" High Duplication Rate", fsize=12), column=1)
        ts.legend.add_face(TextFace("- -", fsize=15, fgcolor="#888888"), column=0)
        ts.legend.add_face(TextFace(" Missing/Uncertain Data", fsize=12), column=1)
        ts.legend_position = 1 # Use integer code 1 for top-right position

        # --- Create Node Styles and Annotations ---
        for node in tree.traverse():
            # Style for all nodes
            ns = NodeStyle()
            ns["size"] = 10
            ns["vt_line_width"] = 2
            ns["hz_line_width"] = 2
            
            # Color branches based on fusion/fission events and stochastic mapping results
            if not node.is_root() and hasattr(node, 'count') and hasattr(node.up, 'count'):
                if node.count > node.up.count:  # Fission
                    ns["hz_line_color"] = "#e41a1c"  # Red
                    ns["hz_line_type"] = 0  # Solid line
                elif node.count < node.up.count:  # Fusion
                    ns["hz_line_color"] = "#377eb8"  # Blue
                    ns["hz_line_type"] = 0  # Solid line
                else:  # Stable
                    ns["hz_line_color"] = "#4daf4a"  # Green
                    ns["hz_line_type"] = 0  # Solid line
                
                # Highlight branches with high duplication rates (from stochastic mapping)
                if hasattr(node, 'expected_duplications') and node.expected_duplications > 0.5:
                    ns["hz_line_color"] = "#ff7f00"  # Orange for high duplication
                    ns["hz_line_width"] = 4  # Thicker line
            else:
                # Default style for branches without count data
                ns["hz_line_color"] = "#888888"  # Gray
                ns["hz_line_type"] = 1  # Dashed line to indicate uncertain/missing data
            
            node.set_style(ns)

            # Add enhanced labels to internal nodes with ML probabilities
            if not node.is_leaf():
                if hasattr(node, 'ml_probability'):
                    label_text = f" n={node.count} (P={node.ml_probability:.2f}) "
                else:
                    label_text = f" n={node.count} "
                
                label = TextFace(label_text, fsize=10, fgcolor="white", bold=True)
                label.margin_top = 2
                label.margin_right = 2
                label.margin_left = 2
                label.margin_bottom = 2
                label.border.width = 1
                label.border.color = "#333333"
                label.background.color = "#666666" # Darker background for contrast
                node.add_face(label, column=0, position="branch-right")
                
                # Add stochastic mapping information if available
                if hasattr(node, 'expected_gains') and hasattr(node, 'expected_losses'):
                    mapping_text = f"G:{node.expected_gains:.1f} L:{node.expected_losses:.1f}"
                    if hasattr(node, 'expected_duplications') and node.expected_duplications > 0.1:
                        mapping_text += f" D:{node.expected_duplications:.1f}"
                    
                    mapping_label = TextFace(mapping_text, fsize=8, fgcolor="#ffcccc")
                    node.add_face(mapping_label, column=0, position="branch-bottom")
            # Leaf names are automatically shown by ts.show_leaf_name = True


        # --- Render and Save Tree ---
        # ETE3 determines format from file extension. DPI is for raster formats.
        render_args = {'tree_style': ts}
        if os.path.splitext(output_file)[1].lower() == '.png':
            render_args['dpi'] = dpi
            
        tree.render(output_file, **render_args)
        
        print(f"* High-quality annotated tree saved to {output_file}")
        print("  - Blue branches indicate fusion events (chromosome number decrease).")
        print("  - Red branches indicate fission events (chromosome number increase).")
        print("  - Green branches indicate stability.")
        print("  - Orange branches indicate high duplication activity.")
        print("  - Gray dashed branches indicate missing or uncertain chromosome data.")
        print("  - Internal nodes show: chromosome number, ML probability, and expected events (G=gains, L=losses, D=duplications).")

    except Exception as e:
        print(f"\n[ERROR] Could not generate the graphical tree. Details: {e}")
        print("  - The analysis and data export parts of the script have completed successfully.")


def analyze_rearrangements(species_a, species_b, mapping_df):
    """
    Analyzes and prints chromosome rearrangements (fissions/fusions) between two species.
    """
    # Check if species exist in the mapping data
    available_species = set(mapping_df['species_A'].unique()) | set(mapping_df['species_B'].unique())
    if species_a not in available_species:
        print(f"\n[WARNING] Species '{species_a}' not found in mapping data. Skipping analysis.")
        return
    if species_b not in available_species:
        print(f"\n[WARNING] Species '{species_b}' not found in mapping data. Skipping analysis.")
        return
    
    # Filter mappings for the two species, considering both directions
    map_ab = mapping_df[((mapping_df['species_A'] == species_a) & (mapping_df['species_B'] == species_b))]
    map_ba = mapping_df[((mapping_df['species_A'] == species_b) & (mapping_df['species_B'] == species_a))]
    
    print(f"\n--- Rearrangement Analysis: {species_a} vs {species_b} ---")

    # Analyze from A to B perspective
    fissions_ab = map_ab[map_ab['bidirectional_mapping_type'].str.contains('Fission', na=False)]
    fusions_ab = map_ab[map_ab['bidirectional_mapping_type'].str.contains('Fusion', na=False)]
    
    if not fissions_ab.empty:
        fission_events = fissions_ab.groupby('chromosome_A')['chromosome_B'].nunique()
        print(f"Fission events (view from {species_a} to {species_b}):")
        for chrom_a, count in fission_events[fission_events > 1].items():
            print(f"  - 1 chromosome in {species_a} ({chrom_a}) splits into {count} in {species_b}.")

    if not fusions_ab.empty:
        fusion_events = fusions_ab.groupby('chromosome_B')['chromosome_A'].nunique()
        print(f"Fusion events (view from {species_a} to {species_b}):")
        for chrom_b, count in fusion_events[fusion_events > 1].items():
            print(f"  - {count} chromosomes in {species_a} fuse into 1 in {species_b} ({chrom_b}).")
            
    if fissions_ab.empty and fusions_ab.empty:
        print("No clear fission or fusion events found in the provided mapping data for this pair.")


def ancestral_reconstruction(tree_file, counts_file, map_file, output_image, output_ancestral_states, output_rearrangements, root_count, layout, show_branch_length, show_support, analyze_pairs, dpi, use_chromevol, optimize_params, model_selection):
    """
    Performs ancestral state reconstruction of chromosome numbers on a phylogenetic tree.
    Enhanced with ChromEvol-inspired maximum likelihood methods.
    """
    # Load the phylogenetic tree
    print(f"--- Loading Data ---")
    print(f"* Tree file: {tree_file}")
    print(f"* Counts file: {counts_file}")
    print(f"* Map file: {map_file}")
    
    # Validate input files exist
    for file_path in [tree_file, counts_file, map_file]:
        if not os.path.exists(file_path):
            print(f"[ERROR] File not found: {file_path}")
            return
    
    tree = Tree(tree_file, format=1)

    # Load the chromosome counts
    counts_df = pd.read_csv(counts_file)
    counts_dict = dict(zip(counts_df['species'], counts_df['chromosome_count']))

    # Attach chromosome counts to the tree leaves
    for leaf in tree:
        if leaf.name in counts_dict:
            leaf.add_features(count=counts_dict[leaf.name])
        else:
            print(f"[WARNING] No chromosome count data found for species: {leaf.name}")

    if use_chromevol:
        # --- ChromEvol-Inspired Maximum Likelihood Analysis ---
        if model_selection:
            model_results = calculate_model_selection_criteria(tree, counts_dict)
        
        # Perform main ChromEvol analysis
        model = perform_chromevol_analysis(tree, counts_dict, 
                                         use_ml=True, 
                                         optimize_params=optimize_params, 
                                         stochastic_mapping=True)
        
        # Update root with constrained value if specified
        root = tree.get_tree_root()
        if root_count is not None:
            print(f"\nConstraining root state to n={root_count} based on user specification.")
            root.add_features(count=root_count)
    else:
        # --- Classical Parsimony-Based Reconstruction ---
        print("\n--- Classical Parsimony-Based Ancestral State Reconstruction ---")
        
        # Phase 1: Post-order traversal (leaves to root)
        for node in tree.traverse("postorder"):
            if not node.is_leaf():
                children_states = [child.states for child in node.children if hasattr(child, 'states')]
                if children_states:
                    intersect_states = set.intersection(*children_states)
                    node.add_features(states=intersect_states if intersect_states else set.union(*children_states))
            else:
                node.add_features(states={node.count} if hasattr(node, 'count') else set())

        # Phase 2: Pre-order traversal
        root = tree.get_tree_root()

        # Constraint based on literature or user specification
        print(f"\nConstraining the root (ancestral) state to n={root_count} based on parameter.")
        root.add_features(count=root_count)
        root.add_features(states={root_count}) # Also update states for consistency

        for node in tree.traverse("preorder"):
            if node.is_root(): continue
            if not node.is_leaf():
                if hasattr(node, 'states') and node.states:
                    parent_states = node.up.states
                    possible_states = list(node.states.intersection(parent_states))
                    node.add_features(count=min(possible_states) if possible_states else min(node.states))
                else:
                    node.add_features(count=node.up.count)
            
            # Only assign auto-generated names to nodes that don't already have names
            if not node.is_leaf() and not node.name:
                child_names = "_".join(sorted([n.name for n in node.get_leaves()]))
                node.name = f"Anc_{len(node.get_leaves())}_{child_names[:30]}"

    print("\nPhylogenetic tree with inferred ancestral chromosome numbers (text version):")
    print(tree.get_ascii(attributes=["name", "count"], show_internal=True))
    
    # Save results
    ancestral_states_data = []
    for i, node in enumerate(tree.traverse()):
        if not node.is_leaf():
            data_row = {
                "node_name": node.name if node.name else f"Internal_Node_{i}", 
                "inferred_chromosome_count": node.count
            }
            
            # Add ML-specific information if available
            if hasattr(node, 'ml_probability'):
                data_row["ml_probability"] = node.ml_probability
            if hasattr(node, 'expected_gains'):
                data_row["expected_gains"] = node.expected_gains
                data_row["expected_losses"] = node.expected_losses
                data_row["expected_duplications"] = node.expected_duplications
            
            ancestral_states_data.append(data_row)
    
    # Save the ancestral states to a CSV file
    ancestral_states_df = pd.DataFrame(ancestral_states_data)
    ancestral_states_df.to_csv(output_ancestral_states, index=False)
    print(f"\n* Ancestral states report saved to {output_ancestral_states}")

    # --- Detailed Analysis of Key Evolutionary Transitions ---
    if analyze_pairs:
        print("\n--- Analysis of Key Chromosomal Rearrangement Events ---")
        mapping_df = pd.read_csv(map_file, sep='\t')
        
        # Process pairs from command line arguments
        if len(analyze_pairs) % 2 != 0:
            print("[WARNING] --analyze_pairs requires an even number of species names (pairs). Ignoring the last unpaired species.")
            analyze_pairs = analyze_pairs[:-1]
        
        for i in range(0, len(analyze_pairs), 2):
            species_a, species_b = analyze_pairs[i], analyze_pairs[i+1]
            analyze_rearrangements(species_a, species_b, mapping_df)
    else:
        print("\n[INFO] No species pairs specified for detailed rearrangement analysis. Use --analyze_pairs to specify pairs.")

    # --- Global Quantification of Rearrangements ---
    quantify_rearrangements_across_tree(tree, output_rearrangements)

    # --- Generate Annotated Tree Image ---
    plot_annotated_tree(tree, output_image, layout, show_branch_length, show_support, dpi)


def quantify_rearrangements_across_tree(tree, output_file):
    """
    Quantifies and reports the number of fusion and fission events across the entire tree.
    """
    print("\n--- Global Analysis of Chromosome Number Evolution ---")
    all_events = []
    for node in tree.traverse():
        if node.is_root() or not hasattr(node, 'count') or not hasattr(node.up, 'count'):
            continue

        parent_count, child_count = node.up.count, node.count
        change = child_count - parent_count
        if change == 0: continue

        event_type = "fission" if change > 0 else "fusion"
        all_events.append({
            "branch": f"Branch from {node.up.name if node.up.name else 'Root'} to {node.name if node.name else 'Unnamed_Node'}",
            "type": event_type,
            "change": abs(change),
            "from_count": parent_count,
            "to_count": child_count
        })

    if not all_events:
        print("No chromosome number changes detected across the tree.")
        return

    if all_events:
        events_df = pd.DataFrame(all_events)
        events_df.sort_values(by=['type', 'change'], ascending=[True, False], inplace=True)
        events_df.to_csv(output_file, index=False)
        print(f"\n* Rearrangement events report saved to {output_file}")


    print("\nDetailed report of chromosomal changes per branch:")
    for event in all_events:
        print(f" - {event['branch']}: {event['type']} of {event['change']} chromosomes (from {event['from_count']} to {event['to_count']})")

    total_fusions = sum(event['change'] for event in all_events if event['type'] == 'fusion')
    total_fissions = sum(event['change'] for event in all_events if event['type'] == 'fission')

    print(f"\nTotal inferred fusion events across the tree: {total_fusions}")
    print(f"Total inferred fission events across the tree: {total_fissions}")

    print("\n--- Overall Evolutionary Trend ---")
    if total_fissions > total_fusions:
        print("The evolutionary history of this group appears to be dominated by chromosome fissions.")
    elif total_fusions > total_fissions:  # Fixed the bug here
        print("The evolutionary history of this group appears to be dominated by chromosome fusions.")
    else:
        print("Chromosome fusions and fissions appear to be relatively balanced across the tree.")


def main():
    """Main function to parse arguments and run the analysis."""
    parser = argparse.ArgumentParser(
        description="Perform ancestral chromosome reconstruction with ChromEvol-inspired methods and visualize the results.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- Input Files ---
    parser.add_argument('--tree', type=str, default='pruned_phylogeny.nwk',
                        help='Path to the phylogenetic tree file in Newick format.')
    parser.add_argument('--counts', type=str, default='chromosome_counts.csv',
                        help='Path to the CSV file with species chromosome counts.')
    parser.add_argument('--map', type=str, default='all.map.tsv',
                        help='Path to the TSV file with chromosome mapping data.')

    # --- Output Files ---
    parser.add_argument('--out_image', type=str, default='annotated_phylogeny.svg',
                        help='Path to save the output phylogeny image.\n'
                             'File extension determines the format (e.g., .svg, .png, .pdf).')
    parser.add_argument('--out_ancestors', type=str, default='ancestral_states_report.csv',
                        help='Path to save the CSV report of inferred ancestral states.')
    parser.add_argument('--out_rearrangements', type=str, default='rearrangement_events_report.csv',
                        help='Path to save the CSV report of rearrangement events.')

    # --- Analysis & Plotting Parameters ---
    parser.add_argument('--root_count', type=int, default=24,
                        help='Constrain the ancestral chromosome number for the root of the tree.')
    parser.add_argument('--layout', type=str, default='circular', choices=['circular', 'rectangular'],
                        help='Tree layout style for visualization.')
    parser.add_argument('--show_branch_length', action='store_true',
                        help='Show branch lengths in the tree visualization.')
    parser.add_argument('--show_support', action='store_true',
                        help='Show branch support values in the tree visualization.')
    parser.add_argument('--analyze_pairs', nargs='*', metavar='SPECIES',
                        help='Specify pairs of species for detailed rearrangement analysis.\n'
                             'Example: --analyze_pairs SpeciesA SpeciesB SpeciesC SpeciesD')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Dots Per Inch (DPI) for raster image formats like PNG.')

    # --- ChromEvol-Inspired Features ---
    parser.add_argument('--use_chromevol', action='store_true',
                        help='Use ChromEvol-inspired maximum likelihood methods instead of parsimony.')
    parser.add_argument('--optimize_params', action='store_true',
                        help='Optimize model parameters using maximum likelihood (requires --use_chromevol).')
    parser.add_argument('--model_selection', action='store_true',
                        help='Perform model selection using AIC/BIC criteria (requires --use_chromevol).')

    args = parser.parse_args()

    ancestral_reconstruction(
        tree_file=args.tree,
        counts_file=args.counts,
        map_file=args.map,
        output_image=args.out_image,
        output_ancestral_states=args.out_ancestors,
        output_rearrangements=args.out_rearrangements,
        root_count=args.root_count,
        layout=args.layout,
        show_branch_length=args.show_branch_length,
        show_support=args.show_support,
        analyze_pairs=args.analyze_pairs,
        dpi=args.dpi,
        use_chromevol=args.use_chromevol,
        optimize_params=args.optimize_params,
        model_selection=args.model_selection
    )


if __name__ == "__main__":
    main()
