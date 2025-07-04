import pandas as pd
from ete3 import Tree
# from ete3.treeview import TreeStyle, NodeStyle, TextFace # Moved into plot_annotated_tree
import argparse
import os
import numpy as np
import scipy.optimize as opt
import scipy.linalg
from scipy.special import logsumexp
# from scipy.stats import poisson # Not used
import itertools # Not used
import warnings
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
warnings.filterwarnings("ignore", category=UserWarning, module="ete3") # Filter ETE3 UserWarnings specifically if needed
warnings.filterwarnings("ignore", category=RuntimeWarning) # Filter general RuntimeWarnings (e.g. from numpy)


# Define model parameter names as a constant for easier maintenance
MODEL_PARAMETER_NAMES = ['gain', 'loss', 'dupl', 'demidupl_gain', 'demidupl_loss', 'base_number_gain', 'base_number_loss', 'linear_rate']

class ChromEvolutionModel:
    """
    ChromEvol-inspired chromosome evolution model with maximum likelihood estimation,
    branch heterogeneity, and parametric modeling of gain, loss, duplication events.
    Includes base number related changes and demi-duplications.
    """

    def __init__(self, max_chromosome_number=50, model_type='full', custom_params=None):
        self.max_chr = max_chromosome_number
        self.model_type = model_type  # 'full', 'simple', 'gain_loss_only'

        # Default model parameters (inspired by ChromEvol)
        # Note: 'demidupl' is split into 'demidupl_gain' and 'demidupl_loss'
        # 'base_number_event_rate' for transitions like i -> i+base_num or i -> i-base_num
        self.params = {
            'gain': 0.1,                # Rate of single chromosome gain (i -> i+1)
            'loss': 0.1,                # Rate of single chromosome loss (i -> i-1)
            'dupl': 0.05,               # Rate of whole genome duplication (i -> 2i)
            'demidupl_gain': 0.01,      # Rate of demi-duplication leading to gain (e.g. i -> 1.5i, then i+k)
            'demidupl_loss': 0.01,      # Rate of demi-duplication leading to loss (e.g. i -> 1.5i, then i-k)
            'base_number': None,        # Base chromosome number (e.g., for x -> x+b or x -> x-b type events)
            'base_number_gain': 0.01,   # Rate of gain by base_number (i -> i + base_number)
            'base_number_loss': 0.01,   # Rate of loss by base_number (i -> i - base_number)
            'linear_rate': 0.01         # Linear chromosome-number dependency factor for gain/loss
        }
        if custom_params: # Allow overriding defaults
            self.params.update(custom_params)

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
            # else:
                # logging.warning(f"Parameter '{key}' not recognized in model.set_parameters().")

    def set_branch_parameters(self, branch_id, **kwargs):
        """Set branch-specific parameters for heterogeneous models"""
        if branch_id not in self.branch_params:
            self.branch_params[branch_id] = self.params.copy()
        for key, value in kwargs.items():
            if key in self.branch_params[branch_id]:
                self.branch_params[branch_id][key] = value

    def build_transition_rate_matrix(self, params=None):
        """
        Build the instantaneous rate matrix Q for chromosome number transitions.
        Based on ChromEvol's approach with gain, loss, duplication, demi-duplication,
        and base number related events.
        """
        current_params = params if params is not None else self.params
        base_number = current_params.get('base_number')

        Q = np.zeros((self.max_chr + 1, self.max_chr + 1))

        for i in range(1, self.max_chr + 1): # Iterate from 1 to max_chr (chromosome number 0 is usually not modeled or has special meaning)
            # Chromosome gain (i -> i+1)
            if i < self.max_chr:
                rate = current_params['gain'] * (1 + current_params['linear_rate'] * i)
                Q[i, i + 1] = rate

            # Chromosome loss (i -> i-1)
            if i > 1:
                rate = current_params['loss'] * (1 + current_params['linear_rate'] * i)
                Q[i, i - 1] = rate

            # Whole genome duplication (i -> 2*i)
            if 2 * i <= self.max_chr and i > 0: # Ensure 2*i is within bounds
                Q[i, 2 * i] += current_params['dupl'] # Use += in case other events also lead to 2*i

            # Demi-duplication (simplified as small fixed gains/losses, could be more complex)
            # Example: i -> 1.5i followed by random loss/gain to integer state
            # For simplicity, we model this as direct transitions to i+1, i+2 (demidupl_gain)
            # and i-1, i-2 (demidupl_loss), if rates are provided.
            # This is a deviation from the original simple /4 distribution if that was the intention.
            # A more ChromEvol-like demi-duplication might involve transitions to i + floor(i/2) or i + ceil(i/2)
            if current_params.get('demidupl_gain', 0) > 0:
                if i + 1 <= self.max_chr:
                    Q[i, i + 1] += current_params['demidupl_gain'] / 2 # Split rate if multiple outcomes
                if i + 2 <= self.max_chr:
                    Q[i, i + 2] += current_params['demidupl_gain'] / 2

            if current_params.get('demidupl_loss', 0) > 0:
                if i - 1 >= 1:
                    Q[i, i - 1] += current_params['demidupl_loss'] / 2
                if i - 2 >= 1:
                    Q[i, i - 2] += current_params['demidupl_loss'] / 2
            
            # Base number related events (e.g., polyploidy followed by chromosome loss leading to x+b or x-b)
            if base_number is not None and base_number > 0:
                # Gain by base_number (i -> i + base_number)
                if i + base_number <= self.max_chr:
                    Q[i, i + base_number] += current_params.get('base_number_gain', 0)

                # Loss by base_number (i -> i - base_number)
                if i - base_number >= 1:
                    Q[i, i - base_number] += current_params.get('base_number_loss', 0)

        # Set diagonal elements: sum of outgoing rates
        for i in range(self.max_chr + 1):
            Q[i, i] = -np.sum(Q[i, :])

        self.Q_matrix = Q
        return Q

    def get_transition_probabilities(self, branch_length, params=None):
        """
        Calculate transition probabilities P(t) = exp(Q*t) for a given branch length.
        """
        current_params = params if params is not None else self.params
        Q = self.build_transition_rate_matrix(current_params)

        if branch_length < 0:
            logging.warning(f"Negative branch length {branch_length} encountered. Using 0.0 instead.")
            branch_length = 0.0
        if branch_length == 0: # If branch length is zero, P is identity matrix
             P = np.eye(self.max_chr + 1)
             return P

        try:
            P = scipy.linalg.expm(Q * branch_length)
            # Ensure probabilities are non-negative and normalized
            P = np.maximum(P, 1e-12) # Increased precision for floor
            P = P / P.sum(axis=1, keepdims=True)
            # Check for NaNs which can arise from expm on ill-conditioned matrices
            if np.isnan(P).any():
                logging.error(f"NaNs in transition probability matrix for branch length {branch_length}. Params: {current_params}")
                # Fallback or error handling needed here. For now, returning identity.
                return np.eye(self.max_chr + 1)
            return P
        except Exception as e:
            logging.error(f"Matrix exponentiation failed for branch length {branch_length}. Error: {e}. Params: {current_params}")
            # Fallback to identity matrix, indicating no change. This is a simplification.
            # A more robust solution might involve trying alternative expm methods or scaling.
            P_fallback = np.eye(self.max_chr + 1)
            return P_fallback

class ChromEvolLikelihoodCalculator:
    """
    Maximum likelihood calculator for chromosome evolution following ChromEvol principles.
    """

    def __init__(self, tree, chromosome_data, model):
        self.tree = tree
        self.data = chromosome_data # species -> count dictionary
        self.model = model
        self.likelihoods = {} # Stores likelihood vectors for each node
        self.root_frequencies = None

    def set_root_frequencies(self, frequencies):
        """Set root state frequencies (prior probabilities for root states)."""
        if frequencies is not None and len(frequencies) == self.model.max_chr + 1:
            self.root_frequencies = frequencies / np.sum(frequencies) # Normalize
        else:
            logging.warning(f"Invalid root frequencies provided or length mismatch. Using uniform distribution. Expected length {self.model.max_chr + 1}, got {len(frequencies) if frequencies is not None else 'None'}")
            self.root_frequencies = np.ones(self.model.max_chr + 1) / (self.model.max_chr + 1)


    def calculate_likelihood_at_node(self, node):
        """
        Calculate likelihood at a given node using Felsenstein's pruning algorithm.
        Likelihood L_i(s) is the probability of observing the data in the clade
        descending from node i, given that node i is in state s.
        """
        node_id = id(node)
        if node_id in self.likelihoods: # Memoization
            return self.likelihoods[node_id]

        likelihood_vector = np.zeros(self.model.max_chr + 1)

        if node.is_leaf():
            observed_count = self.data.get(node.name)
            if observed_count is not None:
                if 0 <= observed_count <= self.model.max_chr:
                    likelihood_vector[observed_count] = 1.0
                else:
                    logging.warning(f"Observed count {observed_count} for leaf {node.name} is outside model range [0, {self.model.max_chr}]. Treating as if count is {self.model.max_chr}.")
                    likelihood_vector[self.model.max_chr] = 1.0 # Or treat as missing data
            else: # Missing data
                logging.info(f"Missing chromosome count data for leaf {node.name}. Treating as equal probability for all states.")
                likelihood_vector.fill(1.0 / (self.model.max_chr + 1)) # Or use specific handling for missing data
            
            self.likelihoods[node_id] = likelihood_vector
            return likelihood_vector

        # Internal node: calculate based on children
        # Initialize with ones for multiplication, or use log-likelihoods
        # Using direct likelihoods here, ensure numerical stability
        current_likelihood_contribution = np.ones(self.model.max_chr + 1)

        for child in node.children:
            child_likelihood_vector = self.calculate_likelihood_at_node(child)
            
            branch_length = child.dist
            if branch_length is None or branch_length < 0: # Handle missing or negative branch lengths
                logging.warning(f"Node {child.name} has invalid branch length {branch_length}. Using default 1e-6.")
                branch_length = 1e-6 # Small positive value to avoid issues with P=I

            # Use branch-specific parameters if available
            P_matrix = self.model.get_transition_probabilities(
                branch_length,
                self.model.branch_params.get(id(child), self.model.params)
            )
            
            # L_parent(s_parent) = sum_{s_child} [ P(s_child | s_parent, t) * L_child(s_child) ]
            # This is effectively (P_matrix @ child_likelihood_vector) for each state of the parent
            conditional_child_likelihood = np.dot(P_matrix, child_likelihood_vector)
            current_likelihood_contribution *= conditional_child_likelihood

            # Normalize to prevent underflow if values get too small
            norm_factor = np.sum(current_likelihood_contribution)
            if norm_factor > 0 and norm_factor < 1e-100 : # Avoid division by zero or near-zero
                 current_likelihood_contribution /= norm_factor
                 # Need to track this scaling factor for the total likelihood if not using log-likelihoods directly.
                 # For simplicity, this example doesn't explicitly track scaling factors for total logL here,
                 # relying on logsumexp later or assuming Felsenstein's original formulation with logs.

        self.likelihoods[node_id] = current_likelihood_contribution
        return current_likelihood_contribution

    def calculate_total_log_likelihood(self):
        """
        Calculate total log-likelihood of the tree given the data and model.
        L = sum_s [ Prior(root_state=s) * L_root(s) ]
        """
        root = self.tree.get_tree_root()
        root_likelihood_vector = self.calculate_likelihood_at_node(root)

        if self.root_frequencies is None: # Ensure root frequencies are set
            self.set_root_frequencies(None) # Sets to uniform if not already set

        # Weighted sum of likelihoods by root state prior probabilities
        # Using logsumexp for numerical stability if root_likelihood_vector contains raw probabilities
        # log_L = logsumexp(np.log(root_likelihood_vector) + np.log(self.root_frequencies))
        # If root_likelihood_vector is already scaled or small, direct sum might be problematic.
        # Assuming calculate_likelihood_at_node returns L_i(s) values that can be directly multiplied.
        
        total_likelihood = np.sum(root_likelihood_vector * self.root_frequencies)

        if total_likelihood <= 0: # Avoid log(0) or log(negative)
            return -np.inf
        return np.log(total_likelihood)

    def ancestral_state_reconstruction(self):
        """
        Perform marginal ancestral state reconstruction (calculates posterior probabilities for each state at each internal node).
        Assigns the most likely state (max posterior probability) to node.ml_count.
        This requires a full upward (for L_i(s)) and downward pass (for posterior).
        """
        if not self.likelihoods: # Ensure up-pass (Felsenstein's) is done
            self.calculate_total_log_likelihood()

        # Downward pass (Felsenstein, 1981 or similar for marginal reconstruction)
        # P(state_k | Data) = L(Data | state_k) * P(state_k) / P(Data)
        # For an internal node u, with parent v and sibling w:
        # Posterior_u(s_u) proportional to Likelihood_u(s_u) * sum_{s_v} [P(s_u|s_v,t_uv) * Posterior_v(s_v) / Likelihood_u_conditional_on_v(s_u)]
        # Simpler: Joint reconstruction often uses the state that maximizes L_root(s) * Prior(s) for the root.
        # For marginal: P(N_i = k | D) = [ (sum_{j} P_ij(t) L_j(k)) * (sum_{m} P_im(t) L_m(k)) * ... ] * pi_k / L_total
        # This is complex. A common approximation for ML state is just the argmax of likelihoods calculated during up-pass for internal nodes.
        # For true marginal probabilities, a more involved calculation is needed.
        # The current code seems to use the up-pass likelihoods directly.

        # Simplified: Assign ML state based on the L_i(s) from the pruning algorithm (up-pass).
        # This is common but is technically the likelihood of the data *below* the node, not the marginal posterior.
        # For proper marginal, one typically does an up-pass then a down-pass.
        
        # Traverse tree (postorder for up-pass is done, pre-order or any order for assigning based on stored likelihoods)
        for node in self.tree.traverse("preorder"): # or "levelorder"
            node_id = id(node)
            if node_id in self.likelihoods:
                node_likelihood_vector = self.likelihoods[node_id]

                # Normalize to get probabilities if they are not already (depends on scaling in up-pass)
                # sum_likelihoods = np.sum(node_likelihood_vector)
                # if sum_likelihoods > 0:
                #     probabilities = node_likelihood_vector / sum_likelihoods
                # else:
                #     probabilities = np.ones_like(node_likelihood_vector) / len(node_likelihood_vector) # Uniform if all zero

                # For now, assume node_likelihood_vector can be used to find max state
                # This is a simplification. True marginal probabilities require more.
                probabilities = node_likelihood_vector / np.sum(node_likelihood_vector) if np.sum(node_likelihood_vector) > 0 else np.ones_like(node_likelihood_vector) / len(node_likelihood_vector)


                most_likely_state = np.argmax(probabilities)
                ml_prob = probabilities[most_likely_state]

                node.add_features(ml_count=most_likely_state)
                node.add_features(ml_probability=ml_prob)
                node.add_features(state_probabilities=probabilities) # Store full distribution

        return self.tree


class ChromEvolOptimizer:
    """
    Parameter optimization for chromosome evolution models using maximum likelihood.
    """
    # Using the constant defined at the module level
    DEFAULT_PARAM_NAMES = [p for p in MODEL_PARAMETER_NAMES if p not in ['base_number']]


    def __init__(self, tree, chromosome_data, model):
        self.tree = tree
        self.data = chromosome_data # species -> count dict
        self.model = model
        # Initialize calculator here, it will be used by objective_function
        self.calculator = ChromEvolLikelihoodCalculator(self.tree, self.data, self.model)


    def objective_function(self, params_vector, param_names_to_optimize):
        """
        Objective function for optimization (negative log-likelihood).
        `param_names_to_optimize` specifies which parameters correspond to `params_vector`.
        """
        try:
            current_model_params = self.model.params.copy()
            update_dict = {}
            for i, param_name in enumerate(param_names_to_optimize):
                # Ensure rates are positive; linear_rate can be small, positive or negative.
                if param_name != 'linear_rate':
                    update_dict[param_name] = max(params_vector[i], 1e-7) # Rates must be positive
                else:
                    update_dict[param_name] = params_vector[i]

            self.model.set_parameters(**update_dict)

            # Reset likelihoods in calculator for new parameter set
            self.calculator.likelihoods = {}
            log_likelihood = self.calculator.calculate_total_log_likelihood()

            # Restore original params if needed, or ensure model is updated only if optimization is successful
            # self.model.set_parameters(**current_model_params) # No, optimizer updates the model
            
            if np.isinf(log_likelihood) or np.isnan(log_likelihood):
                return 1e9 # Large penalty for invalid likelihood

            return -log_likelihood  # Return negative for minimization
        except Exception as e:
            logging.error(f"Error in objective function: {e} with params {params_vector}")
            return 1e10 # Large penalty for any other error

    def optimize_parameters(self, method='L-BFGS-B', maxiter=1000, params_to_optimize=None):
        """
        Optimize model parameters using maximum likelihood.
        `params_to_optimize`: A list of parameter names to optimize. If None, uses default.
        """
        if params_to_optimize is None:
            params_to_optimize = self.DEFAULT_PARAM_NAMES
        
        # Filter out params not in the model or non-numeric (like base_number if it's fixed or None)
        valid_params_to_optimize = [p for p in params_to_optimize if p in self.model.params and isinstance(self.model.params[p], (int, float))]

        initial_params_vector = [self.model.params[p_name] for p_name in valid_params_to_optimize]

        # Parameter bounds (rates generally > 0, linear_rate can be negative)
        # Bounds should correspond to `valid_params_to_optimize`
        bounds = []
        for p_name in valid_params_to_optimize:
            if p_name == 'linear_rate':
                bounds.append((-1.0, 1.0)) # Example: linear rate bounded
            elif p_name == 'base_number': # Should not be optimized directly if it's a structural param
                continue # Or handle as fixed
            else: # Rates
                bounds.append((1e-7, 10.0)) # Rates positive and bounded

        logging.info(f"Optimizing model parameters using {method} for: {valid_params_to_optimize}")
        logging.info(f"Initial parameters: {dict(zip(valid_params_to_optimize, initial_params_vector))}")
        
        if not initial_params_vector:
            logging.warning("No parameters to optimize.")
            return None

        # Optimize
        result = opt.minimize(
            self.objective_function,
            np.array(initial_params_vector), # Ensure it's a numpy array
            args=(valid_params_to_optimize,),
            method=method,
            bounds=bounds,
            options={'maxiter': maxiter, 'disp': False} # disp=True for more verbose output
        )

        if result.success:
            optimized_values = result.x
            final_params_dict = dict(zip(valid_params_to_optimize, optimized_values))
            self.model.set_parameters(**final_params_dict) # Update model with optimized parameters
            
            logging.info("Optimization successful!")
            log_params = {k: f"{v:.4f}" for k, v in final_params_dict.items()}
            logging.info(f"Optimized parameters: {log_params}")
            logging.info(f"Final log-likelihood: {-result.fun:.4f}")
            return result
        else:
            logging.error(f"Optimization failed: {result.message}")
            return None

class StochasticMapper:
    """
    Stochastic mapping implementation for chromosome evolution events.
    Simulates actual evolutionary histories (event timings and types) along branches
    consistent with the model and inferred ancestral states.
    """

    def __init__(self, tree, model, num_mappings=100):
        self.tree = tree # Tree should have ancestral states (ml_count) inferred
        self.model = model
        self.num_mappings = num_mappings
        # self.event_counts = {} # This will store summary of events per branch

    def _classify_event(self, from_state, to_state, params):
        """Helper to classify event type based on state change and model parameters."""
        if to_state == from_state + 1: return 'gain'
        if to_state == from_state - 1: return 'loss'
        if to_state == 2 * from_state and from_state > 0: return 'duplication'

        base_num = params.get('base_number')
        if base_num and base_num > 0:
            if to_state == from_state + base_num: return 'base_number_gain'
            if to_state == from_state - base_num: return 'base_number_loss'

        # More complex demi-duplication events might be harder to classify directly
        # if they result in various i+k states.
        # For now, 'other' covers these or unclassified transitions.
        # if abs(to_state - 1.5 * from_state) < 2 : return 'demi-duplication_related' # Example heuristic

        # Check for demi-duplication outcomes if they are distinct enough
        # This part is tricky as demidupl rates affect multiple transitions (e.g. i+1, i+2)
        # For simplicity, changes not fitting above are 'other'
        delta = to_state - from_state
        if delta > 1 and delta < from_state : return 'demi_gain_like' # Heuristic
        if delta < -1 and abs(delta) < from_state/2 : return 'demi_loss_like' # Heuristic

        return 'other'


    def simulate_events_on_branch(self, start_state, branch_length, params):
        """
        Simulate chromosome evolution events along a single branch using Gillespie-like approach.
        The simulation proceeds for the given branch_length, not necessarily until an end_state is reached.
        """
        events = []
        current_state = start_state
        current_time = 0.0

        Q = self.model.build_transition_rate_matrix(params)

        while current_time < branch_length:
            if current_state < 0 or current_state > self.model.max_chr:
                logging.warning(f"Simulate_events: current_state {current_state} out of bounds. Stopping simulation on branch.")
                break

            total_rate_out = -Q[current_state, current_state]
            if total_rate_out <= 1e-9: # Effectively zero rate out
                current_time = branch_length # No more events will occur, advance time to end
                break

            # Time to next event
            time_to_event = np.random.exponential(1.0 / total_rate_out)

            if current_time + time_to_event >= branch_length:
                current_time = branch_length # Event occurs after branch ends
                break
            
            current_time += time_to_event

            # Choose next state based on transition rates from current_state
            # Probabilities of transitioning to state j are Q[i,j] / total_rate_out
            possible_next_states_rates = Q[current_state, :].copy()
            possible_next_states_rates[current_state] = 0 # Do not transition to self here
            
            # Ensure rates are non-negative before normalization
            possible_next_states_rates = np.maximum(possible_next_states_rates, 0)

            if np.sum(possible_next_states_rates) <= 1e-9: # No valid transitions out
                 current_time = branch_length # Should have been caught by total_rate_out check
                 break

            probabilities = possible_next_states_rates / np.sum(possible_next_states_rates)
            next_state = np.random.choice(np.arange(self.model.max_chr + 1), p=probabilities)

            event_type = self._classify_event(current_state, next_state, params)
            events.append({
                'time_on_branch': current_time,
                'from_state': current_state,
                'to_state': next_state,
                'type': event_type
            })
            current_state = next_state
        
        return events, current_state # Return final state on branch as well

    def perform_stochastic_mapping(self):
        """
        Perform stochastic mapping across the entire tree for `self.num_mappings` replicates.
        Requires ancestral states (node.ml_count) to be populated on the tree.
        """
        all_branch_mappings = {id(node): [] for node in self.tree.traverse() if not node.is_root()}

        logging.info(f"Performing stochastic mapping with {self.num_mappings} replicates...")

        for i_mapping in range(self.num_mappings):
            if (i_mapping + 1) % (self.num_mappings // 10 or 1) == 0: # Log progress
                logging.info(f"  Stochastic mapping replicate {i_mapping + 1}/{self.num_mappings}")

            # For each replicate, simulate history down the tree from the root
            # Root state for simulation can be drawn from posterior or fixed to ML state
            root = self.tree.get_tree_root()
            
            # Initialize states for this replicate's traversal
            # node_states_this_replicate = {id(root): root.ml_count} # Start with ML state at root
            
            # Option: Draw root state from its posterior probability distribution if available
            if hasattr(root, 'state_probabilities'):
                current_root_state = np.random.choice(np.arange(self.model.max_chr + 1), p=root.state_probabilities)
            else:
                current_root_state = root.ml_count # Fallback to ML state
            node_states_this_replicate = {id(root): current_root_state}


            for node in self.tree.traverse("preorder"): # Traverse from root downwards
                if node.is_root():
                    continue

                parent_id = id(node.up)
                start_state_for_branch = node_states_this_replicate.get(parent_id)

                if start_state_for_branch is None: # Should not happen in preorder traversal if root is handled
                    logging.error(f"Parent state not found for node {node.name} in stochastic mapping replicate {i_mapping+1}. Skipping branch.")
                    continue

                branch_length = node.dist
                if branch_length is None or branch_length <= 0: # Skip zero/negative length branches for simulation
                    node_states_this_replicate[id(node)] = start_state_for_branch # State persists
                    all_branch_mappings[id(node)].append([]) # No events on this zero-length branch
                    continue

                # Branch specific params or global model params
                branch_params = self.model.branch_params.get(id(node), self.model.params)

                branch_events, end_state_on_branch = self.simulate_events_on_branch(
                    start_state_for_branch, branch_length, branch_params
                )

                node_states_this_replicate[id(node)] = end_state_on_branch # Store simulated end state for children
                all_branch_mappings[id(node)].append(branch_events) # Store events for this replicate

        # Summarize event counts from all_branch_mappings
        self.event_counts_summary = self._summarize_event_counts(all_branch_mappings)
        
        # Annotate tree with mean event counts
        for node in self.tree.traverse():
            if not node.is_root():
                branch_id = id(node)
                if branch_id in self.event_counts_summary:
                    summary = self.event_counts_summary[branch_id]
                    for event_type, metrics in summary.items():
                        node.add_features(**{f"expected_{event_type}": metrics['mean']})

        return self.event_counts_summary


    def _summarize_event_counts(self, all_branch_mappings_across_replicates):
        """
        Summarize event counts for each branch across all stochastic mappings.
        `all_branch_mappings_across_replicates`: dict {branch_id: [[events_rep1], [events_rep2], ...]}
        """
        summary_per_branch = {}

        # Define all event types we want to track, including new ones
        # This should align with _classify_event outputs
        event_types_to_track = ['gain', 'loss', 'duplication',
                                'base_number_gain', 'base_number_loss',
                                'demi_gain_like', 'demi_loss_like', 'other']


        for branch_id, mappings_for_branch in all_branch_mappings_across_replicates.items():
            # mappings_for_branch is a list of lists of event dicts (one list per replicate)

            branch_summary = {etype: {'counts_all_reps': []} for etype in event_types_to_track}

            for single_replicate_events in mappings_for_branch:
                # Count events of each type for this replicate
                counts_this_replicate = {etype: 0 for etype in event_types_to_track}
                for event in single_replicate_events:
                    etype = event['type']
                    if etype in counts_this_replicate:
                        counts_this_replicate[etype] += 1
                    else: # Should not happen if event_types_to_track is comprehensive
                        counts_this_replicate.setdefault('unknown', 0)
                        counts_this_replicate['unknown'] +=1
                
                # Store these counts
                for etype in event_types_to_track:
                    branch_summary[etype]['counts_all_reps'].append(counts_this_replicate[etype])

            # Calculate mean, std, CI for each event type on this branch
            final_branch_metrics = {}
            for etype, data in branch_summary.items():
                counts = data['counts_all_reps']
                if counts: # Ensure there are counts to process
                    final_branch_metrics[etype] = {
                        'mean': np.mean(counts),
                        'std': np.std(counts),
                        'ci_lower': np.percentile(counts, 2.5),
                        'ci_upper': np.percentile(counts, 97.5),
                        'total_occurrences': np.sum(counts) # Total across all replicates
                    }
                else: # No events of this type or no mappings for this branch
                     final_branch_metrics[etype] = {'mean': 0, 'std': 0, 'ci_lower': 0, 'ci_upper': 0, 'total_occurrences':0}

            summary_per_branch[branch_id] = final_branch_metrics

        return summary_per_branch

# Enhanced analysis functions

def perform_chromevol_analysis(tree, counts_dict, max_chr_num, root_prior_state=None, use_ml=True, optimize_params=True, stochastic_mapping=True, model_params_input=None):
    """
    Perform ChromEvol-inspired maximum likelihood analysis of chromosome evolution.
    `max_chr_num`: Maximum chromosome number to model.
    `root_prior_state`: Integer for a fixed root state, or None for uniform/estimated.
    `model_params_input`: Dictionary of initial model parameters.
    """
    logging.info("\n--- ChromEvol-Inspired Maximum Likelihood Analysis ---")

    # Initialize model
    model = ChromEvolutionModel(max_chromosome_number=max_chr_num, model_type='full', custom_params=model_params_input)
    
    # Attach data to tree leaves
    # Ensure all leaves have a 'count' attribute, even if None (for missing data handling in LikelihoodCalculator)
    max_obs_count = 0
    for leaf in tree:
        if leaf.name in counts_dict:
            count = counts_dict[leaf.name]
            leaf.add_features(count=count)
            if count is not None:
                 max_obs_count = max(max_obs_count, count)
        else:
            leaf.add_features(count=None) # Explicitly mark missing data
            logging.warning(f"No chromosome count for species {leaf.name} in input file.")

    if max_obs_count > model.max_chr:
        logging.warning(f"Max observed chromosome count ({max_obs_count}) > model.max_chr ({model.max_chr}). Consider increasing max_chr_num.")


    if use_ml:
        calculator = ChromEvolLikelihoodCalculator(tree, counts_dict, model)
        
        # Set root frequencies
        if root_prior_state is not None and 0 <= root_prior_state <= model.max_chr:
            logging.info(f"Setting strong prior for root state = {root_prior_state}")
            root_freq = np.zeros(model.max_chr + 1)
            root_freq[root_prior_state] = 1.0
        else:
            logging.info("Using uniform prior for root state frequencies.")
            root_freq = np.ones(model.max_chr + 1) / (model.max_chr + 1)
        calculator.set_root_frequencies(root_freq)

        if optimize_params:
            optimizer = ChromEvolOptimizer(tree, counts_dict, model) # Pass calculator or re-init
            optimizer.calculator = calculator # Ensure optimizer uses the same calculator instance
            result = optimizer.optimize_parameters() # Uses default list of params to optimize
            if not result or not result.success:
                logging.error("Parameter optimization failed. Continuing with initial/default parameters.")
        
        # Calculate final likelihood with current (possibly optimized) parameters
        # Need to reset likelihoods if parameters changed and calculator was not part of optimizer's direct update loop
        calculator.likelihoods = {} # Reset before final calculation
        final_log_likelihood = calculator.calculate_total_log_likelihood()
        logging.info(f"Final log-likelihood with current parameters: {final_log_likelihood:.4f}")

        # Perform ancestral reconstruction
        calculator.ancestral_state_reconstruction() # This populates ml_count, ml_probability on nodes

        # Copy ML results to a standard 'count' attribute for visualization consistency if needed by plotting
        for node in tree.traverse():
            if hasattr(node, 'ml_count'): # Check if ChromEvol path was taken
                node.add_features(count=node.ml_count) # Overwrite/set general 'count'
        
        if stochastic_mapping:
            mapper = StochasticMapper(tree, model, num_mappings=100) # Tree now has ml_counts
            event_summary = mapper.perform_stochastic_mapping() # Populates expected_event features
            
            logging.info("\n--- Expected Number of Events per Branch (Stochastic Mapping Summary) ---")
            for node in tree.traverse():
                if not node.is_root() and hasattr(node, 'expected_gain'): # Check one example attribute
                    node_id = id(node)
                    branch_events_summary = mapper.event_counts_summary.get(node_id, {})
                    log_msg_parts = []
                    for etype, metrics in branch_events_summary.items():
                        if metrics['mean'] > 0.01 : # Log if mean is somewhat significant
                            log_msg_parts.append(f"{etype}={metrics['mean']:.2f}")
                    if log_msg_parts:
                         node_name_str = node.name if node.name else f"Branch_to_{node.get_leaves()[0].name[:10]}"
                         logging.info(f"  {node_name_str}: {', '.join(log_msg_parts)}")
    else:
        # This implies parsimony or other non-ML method, not fully implemented here with 'use_ml=False'
        logging.info("use_ml=False, skipping ML-based analysis and stochastic mapping.")
        # Fallback to a simple parsimony if that's intended, or error if ML is required
        # The original script had a parsimony path; this refactor focuses on ML.
        # For now, if not use_ml, no ChromEvol specific analysis is done.
        pass
    
    return model, tree # Return model and the (potentially) annotated tree

def calculate_model_selection_criteria(tree, counts_dict, max_chr_num, root_prior_state=None, model_configurations=None):
    """
    Calculate AIC/BIC for different model configurations.
    `model_configurations`: A list of dicts, each specifying a model setup (e.g., params to fix/optimize).
    Example: [{'name': 'simple_model', 'params_to_optimize': ['gain', 'loss']},
              {'name': 'full_model', 'params_to_optimize': ChromEvolOptimizer.DEFAULT_PARAM_NAMES}]
    """
    logging.info("\n--- Model Selection based on AIC/BIC ---")
    
    results = {}
    
    if model_configurations is None: # Default: compare a simple vs. a fuller model
        model_configurations = [
            {'name': 'GainLossOnly', 'initial_params': {'gain':0.1, 'loss':0.1, 'dupl':0, 'demidupl_gain':0, 'demidupl_loss':0, 'base_number_gain':0, 'base_number_loss':0, 'linear_rate':0},
             'params_to_optimize': ['gain', 'loss']},
            {'name': 'FullModelDefaultRates', 'initial_params': None, # Uses class defaults
             'params_to_optimize': ChromEvolOptimizer.DEFAULT_PARAM_NAMES}, # Optimize all standard rates
        ]

    for config in model_configurations:
        model_name = config['name']
        initial_params_for_model = config.get('initial_params') # Can be None to use class defaults
        params_to_optimize_for_model = config['params_to_optimize']
        
        logging.info(f"Evaluating model: {model_name} with params to optimize: {params_to_optimize_for_model}")

        # Initialize model with specific configuration if provided
        model = ChromEvolutionModel(max_chromosome_number=max_chr_num, custom_params=initial_params_for_model)
        
        # Setup calculator
        calculator = ChromEvolLikelihoodCalculator(tree, counts_dict, model)
        if root_prior_state is not None and 0 <= root_prior_state <= model.max_chr:
            root_freq = np.zeros(model.max_chr + 1); root_freq[root_prior_state] = 1.0
        else:
            root_freq = np.ones(model.max_chr + 1) / (model.max_chr + 1)
        calculator.set_root_frequencies(root_freq)

        # Optimizer
        optimizer = ChromEvolOptimizer(tree, counts_dict, model)
        optimizer.calculator = calculator # Share calculator instance

        opt_result = optimizer.optimize_parameters(params_to_optimize=params_to_optimize_for_model, maxiter=500)
        
        if opt_result and opt_result.success:
            log_likelihood = -opt_result.fun
            # Number of parameters optimized is len(params_to_optimize_for_model)
            # Ensure base_number is not counted if it's fixed/structural
            n_params_optimized = len([p for p in params_to_optimize_for_model if p != 'base_number'])
            
            n_datapoints = len([leaf for leaf in tree if leaf.name in counts_dict and counts_dict[leaf.name] is not None]) # Number of species with data

            AIC = 2 * n_params_optimized - 2 * log_likelihood
            BIC = np.log(n_datapoints) * n_params_optimized - 2 * log_likelihood if n_datapoints > 0 else AIC # Avoid log(0)

            results[model_name] = {
                'log_likelihood': log_likelihood,
                'AIC': AIC,
                'BIC': BIC,
                'n_params': n_params_optimized,
                'optimized_parameters': {k:v for k,v in model.params.items() if k in params_to_optimize_for_model}
            }
            logging.info(f"  Model {model_name}: LogL={log_likelihood:.2f}, AIC={AIC:.2f}, BIC={BIC:.2f}, N_params={n_params_optimized}")
        else:
            logging.warning(f"  Optimization failed for model: {model_name}")
            results[model_name] = {'log_likelihood': -np.inf, 'AIC': np.inf, 'BIC': np.inf, 'n_params': 0, 'optimized_parameters': {}}

    if results:
        try:
            best_aic_model = min(results.keys(), key=lambda k: results[k]['AIC'])
            best_bic_model = min(results.keys(), key=lambda k: results[k]['BIC'])
            logging.info(f"\nBest model by AIC: {best_aic_model} (AIC: {results[best_aic_model]['AIC']:.2f})")
            logging.info(f"Best model by BIC: {best_bic_model} (BIC: {results[best_bic_model]['BIC']:.2f})")
        except ValueError: # Handles cases where all models might have failed and results are empty or contain infs
            logging.error("Could not determine best model due to invalid AIC/BIC values.")

    return results


def plot_annotated_tree(tree, output_file, layout, show_branch_length, show_support, dpi, title="Ancestral Chromosome Evolution"):
    """
    Generates and saves an image of the phylogenetic tree with annotations.
    Tree should have 'count' for chromosome numbers, 'ml_probability', and 'expected_...' for events.
    """
    logging.info("\n--- Generating Annotated Phylogenetic Tree Image ---")
    try:
        from ete3.treeview import TreeStyle, NodeStyle, TextFace # Import here
    except ImportError as e:
        logging.warning(f"Failed to import ete3.treeview components, likely due to missing GUI dependencies: {e}")
        logging.warning("Skipping tree image generation. Other analysis results should still be available.")
        return # Skip plotting

    try:
        ts = TreeStyle()
        ts.mode = "c" if layout == "circular" else "r"
        ts.show_leaf_name = True
        ts.show_branch_length = show_branch_length
        ts.show_branch_support = show_support
        ts.scale = 120 # Adjust as needed
        ts.title.add_face(TextFace(title, fsize=15, fgcolor="black", bold=True), column=0)

        # Basic legend (can be expanded)
        # Colors used: Fusion=Blue (#377eb8), Fission=Red (#e41a1c), Stable=Green (#4daf4a), HighDupl=Orange (#ff7f00)
        legend_elements = [
            ("■", "#377eb8", "Fusion (n decreases)"),
            ("■", "#e41a1c", "Fission (n increases)"),
            ("■", "#4daf4a", "Stable Lineage"),
            ("■", "#ff7f00", "High Duplication Rate"),
            ("- -", "#888888", "Uncertain/Missing Data")
        ]
        for symbol, color, text in legend_elements:
            ts.legend.add_face(TextFace(symbol, fsize=15, fgcolor=color), column=0)
            ts.legend.add_face(TextFace(f" {text}", fsize=12), column=1)
        ts.legend_position = 1 # Top-right

        for node in tree.traverse():
            ns = NodeStyle()
            ns["size"] = 8 if node.is_leaf() else 10
            ns["vt_line_width"] = 2
            ns["hz_line_width"] = 2
            ns["shape"] = "sphere" if not node.is_leaf() else "circle"


            # Branch coloring based on chromosome number change (fission/fusion)
            # And highlighting for high duplication if data is available
            if not node.is_root() and hasattr(node, 'count') and node.count is not None and \
               hasattr(node.up, 'count') and node.up.count is not None:
                
                parent_count = node.up.count
                child_count = node.count

                if child_count > parent_count: ns["hz_line_color"] = "#e41a1c" # Fission (Red)
                elif child_count < parent_count: ns["hz_line_color"] = "#377eb8" # Fusion (Blue)
                else: ns["hz_line_color"] = "#4daf4a" # Stable (Green)

                # Highlight for high duplication (e.g., expected_duplication > threshold)
                # This threshold (0.5) is arbitrary, adjust based on typical values from stochastic mapping
                if hasattr(node, 'expected_duplication') and node.expected_duplication > 0.5:
                    ns["hz_line_color"] = "#ff7f00" # Orange
                    ns["hz_line_width"] = 4 # Thicker line
            else: # Root branch or missing data
                ns["hz_line_color"] = "#888888" # Gray
                if not node.is_leaf() : ns["hz_line_type"] = 1 # Dashed for internal uncertain branches


            node.set_style(ns)

            # Annotations for internal nodes: chromosome count, ML probability, event expectations
            if not node.is_leaf():
                label_parts = []
                if hasattr(node, 'count') and node.count is not None:
                    label_parts.append(f"n={node.count}")
                if hasattr(node, 'ml_probability') and node.ml_probability is not None:
                     label_parts.append(f"(P={node.ml_probability:.2f})")
                
                if label_parts:
                    node_label_text = " ".join(label_parts)
                    label_face = TextFace(node_label_text, fsize=9, fgcolor="white", bold=True)
                    label_face.background.color = "#555555"
                    label_face.margin_top = 2; label_face.margin_bottom = 2
                    label_face.margin_left = 2; label_face.margin_right = 2
                    label_face.border.width = 1; label_face.border.color = "#333333"
                    node.add_face(label_face, column=0, position="branch-right")

                # Stochastic mapping annotations (example: Gains, Losses, Duplications)
                event_texts = []
                if hasattr(node, 'expected_gain') and node.expected_gain > 0.05: event_texts.append(f"G:{node.expected_gain:.1f}")
                if hasattr(node, 'expected_loss') and node.expected_loss > 0.05: event_texts.append(f"L:{node.expected_loss:.1f}")
                if hasattr(node, 'expected_duplication') and node.expected_duplication > 0.05: event_texts.append(f"D:{node.expected_duplication:.1f}")
                
                if event_texts:
                    event_label_text = " ".join(event_texts)
                    event_face = TextFace(event_label_text, fsize=7, fgcolor="#dddddd") # Light text for dark branches
                    node.add_face(event_face, column=0, position="branch-bottom")
            else: # Leaf nodes
                # ns["fgcolor"] = "black" # Leaf name color if needed
                pass

        render_args = {'tree_style': ts}
        if os.path.splitext(output_file)[1].lower() == '.png':
            render_args['dpi'] = dpi
        tree.render(output_file, **render_args)
        logging.info(f"* High-quality annotated tree saved to {output_file}")

    except ImportError: # Already caught above, but good to have a specific catch if other imports were here
        # This block might be redundant if the first try-except for imports handles all relevant ImportErrors
        logging.warning("Skipping tree image generation due to import error (already logged).")
    except Exception as e: # Catch other errors during tree rendering (e.g., display errors, font issues)
        logging.warning(f"Could not generate the graphical tree, possibly due to missing display or rendering libraries. Details: {e}")
        # logging.debug("Detailed error for tree generation:", exc_info=True) # Optional: more detail for debugging
        logging.info("  - Skipping tree image generation. Other analysis results should still be available.")


def analyze_rearrangements(species_a, species_b, mapping_df):
    """
    Analyzes and prints chromosome rearrangements (fissions/fusions) between two species
    based on a synteny/collinearity mapping dataframe.
    """
    available_species = set(mapping_df['species_A'].unique()) | set(mapping_df['species_B'].unique())
    if species_a not in available_species:
        logging.warning(f"Species '{species_a}' not found in mapping data. Skipping analysis for this pair.")
        return
    if species_b not in available_species:
        logging.warning(f"Species '{species_b}' not found in mapping data. Skipping analysis for this pair.")
        return

    # Filter mappings for the two species, considering both A-B and B-A mappings if columns are fixed
    # Assuming 'bidirectional_mapping_type' column exists and indicates Fission/Fusion from perspective of species_A
    map_ab = mapping_df[((mapping_df['species_A'] == species_a) & (mapping_df['species_B'] == species_b))]
    
    logging.info(f"\n--- Rearrangement Analysis: {species_a} (A) vs {species_b} (B) ---")

    # Fissions: 1 chr in A maps to N chrs in B
    fissions_from_a = map_ab[map_ab['bidirectional_mapping_type'].str.contains('Fission', na=False)]
    if not fissions_from_a.empty:
        # Count how many B chromosomes each A chromosome maps to
        fission_events = fissions_from_a.groupby('chromosome_A')['chromosome_B'].nunique()
        logging.info(f"Fission events (1 chr in {species_a} -> N chrs in {species_b}):")
        for chrom_a, num_b_chrs in fission_events[fission_events > 1].items(): # N > 1 indicates fission
            logging.info(f"  - {species_a} chr {chrom_a} splits into {num_b_chrs} chromosomes in {species_b}.")

    # Fusions: N chrs in A map to 1 chr in B
    fusions_to_b = map_ab[map_ab['bidirectional_mapping_type'].str.contains('Fusion', na=False)]
    if not fusions_to_b.empty:
        # Count how many A chromosomes map to each B chromosome
        fusion_events = fusions_to_b.groupby('chromosome_B')['chromosome_A'].nunique()
        logging.info(f"Fusion events (N chrs in {species_a} -> 1 chr in {species_b}):")
        for chrom_b, num_a_chrs in fusion_events[fusion_events > 1].items(): # N > 1 indicates fusion
            logging.info(f"  - {num_a_chrs} chromosomes in {species_a} fuse into {species_b} chr {chrom_b}.")
            
    if fissions_from_a.empty and fusions_to_b.empty:
        logging.info("No clear fission or fusion events found in the provided mapping data for this pair and direction.")


def ancestral_reconstruction(tree_file, counts_file, map_file,
                             output_image, output_ancestral_states, output_rearrangements,
                             root_constraint_count, layout, show_branch_length, show_support,
                             analyze_pairs_list, dpi,
                             use_chromevol_method, optimize_model_params, perform_model_selection,
                             max_chromosome_val, initial_model_params_dict): # New params
    """
    Main pipeline for ancestral state reconstruction of chromosome numbers.
    Enhanced with ChromEvol-inspired maximum likelihood methods.
    """
    logging.info(f"--- Loading Data ---")
    logging.info(f"* Tree file: {tree_file}")
    logging.info(f"* Counts file: {counts_file}")
    if map_file: logging.info(f"* Map file: {map_file}")

    # Validate essential input files
    if not os.path.exists(tree_file): logging.error(f"Tree file not found: {tree_file}"); return
    if not os.path.exists(counts_file): logging.error(f"Counts file not found: {counts_file}"); return
    if map_file and not os.path.exists(map_file): logging.warning(f"Map file specified but not found: {map_file}")

    tree = Tree(tree_file, format=1) # Assuming Newick format 1 as default for ETE3
    counts_df = pd.read_csv(counts_file)
    counts_dict = dict(zip(counts_df['species'], counts_df['chromosome_count']))

    # Determine max_chr dynamically if not provided or if current one is too low
    observed_max_count = counts_df['chromosome_count'].max()
    if max_chromosome_val is None:
        max_chromosome_val = observed_max_count + 10 # Add a buffer
        logging.info(f"Dynamically set max_chromosome_val to {max_chromosome_val} based on input data.")
    elif observed_max_count > max_chromosome_val :
        logging.warning(f"Max observed chromosome count ({observed_max_count}) is greater than specified max_chromosome_val ({max_chromosome_val}). Consider increasing it.")


    # Attach initial counts to leaves (will be used by both ML and Parsimony paths)
    for leaf in tree:
        if leaf.name in counts_dict:
            leaf.add_features(count=counts_dict[leaf.name])
        else:
            logging.warning(f"No chromosome count data found for species: {leaf.name}. Will be treated as missing data.")
            leaf.add_features(count=None)


    if use_chromevol_method:
        logging.info("\n--- Using ChromEvol-Inspired Maximum Likelihood Method ---")
        
        if perform_model_selection:
            # This would typically run perform_chromevol_analysis multiple times with different model settings
            # For now, assume it returns the best model or parameters to use
            model_sel_results = calculate_model_selection_criteria(tree, counts_dict, max_chromosome_val, root_constraint_count)
            # Logic to pick the best model and its parameters from model_sel_results would go here.
            # For simplicity, we'll just proceed with the standard analysis after printing selection results.
            # Ideally, the chosen model's parameters would be passed to perform_chromevol_analysis.
            logging.info(f"Model selection results: {model_sel_results}")


        # Main ML analysis
        # Pass root_constraint_count as root_prior_state
        # Pass initial_model_params_dict for model initialization
        final_model, annotated_tree = perform_chromevol_analysis(
            tree, counts_dict, max_chromosome_val,
            root_prior_state=root_constraint_count,
            use_ml=True, # This is implicit in use_chromevol_method
            optimize_params=optimize_model_params,
            stochastic_mapping=True, # Assuming always true if ChromEvol is used, can be a param
            model_params_input=initial_model_params_dict
        )
        tree = annotated_tree # Update tree with annotations from ML analysis

    else:
        # --- Classical Parsimony-Based Reconstruction ---
        logging.info("\n--- Using Classical Parsimony-Based Ancestral State Reconstruction ---")
        # (This part of the original code can be kept if parsimony is still a desired option)
        # Phase 1: Post-order traversal (Sankoff-like, finding sets of possible states)
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                if hasattr(node, 'count') and node.count is not None:
                    node.add_features(states={node.count})
                else:
                    node.add_features(states=set(range(1, max_chromosome_val + 1))) # All states possible if missing
            else: # Internal node
                children_states_list = [child.states for child in node.children if hasattr(child, 'states')]
                if children_states_list:
                    # Intersection of children's states (most parsimonious)
                    intersect_states = set.intersection(*children_states_list)
                    if intersect_states:
                        node.add_features(states=intersect_states)
                    else: # If intersection is empty, take union (less parsimonious, implies change)
                        node.add_features(states=set.union(*children_states_list))
                else: # No children with states (should not happen for internal nodes in postorder)
                     node.add_features(states=set(range(1, max_chromosome_val + 1)))

        # Phase 2: Pre-order traversal (assigning a single state)
        root = tree.get_tree_root()
        if root_constraint_count is not None and root_constraint_count in root.states:
            root.add_features(count=root_constraint_count)
        elif root.states: # Pick a state (e.g., min, or one minimizing changes to children)
            root.add_features(count=min(root.states)) # Simple choice: min state
        else: # Should not happen
            root.add_features(count=max_chromosome_val // 2)
        root.states = {root.count} # Finalize root state

        for node in tree.traverse("preorder"):
            if node.is_root(): continue # Root already handled

            if not hasattr(node, 'count'): # If count not assigned yet
                parent_state = node.up.count
                if parent_state in node.states: # Prefer parent state if possible (no change)
                    node.add_features(count=parent_state)
                elif node.states: # Choose state from its possible set (e.g., min, or closest to parent)
                    # Choose state in node.states closest to parent_state
                    node.add_features(count=min(node.states, key=lambda s: abs(s - parent_state)))
                else: # Fallback
                    node.add_features(count=parent_state)
            node.states = {node.count} # Finalize node state
            
            if not node.is_leaf() and not node.name: # Assign generic name if missing
                child_leaves = node.get_leaves()
                child_names_str = "_".join(sorted([n.name for n in child_leaves])[:3]) # Max 3 names
                node.name = f"Anc_{len(child_leaves)}_{child_names_str[:20]}" # Ensure name not too long


    logging.info("\nPhylogenetic tree with inferred ancestral chromosome numbers (text version):")
    try:
        logging.info(tree.get_ascii(attributes=["name", "count", "ml_probability"], show_internal=True))
    except Exception as e:
        logging.warning(f"Could not print tree ASCII representation: {e}")

    # --- Save Results ---
    ancestral_states_data_list = []
    for i, node in enumerate(tree.traverse()):
        # Save data for internal nodes and optionally for leaves if they were processed
        # if not node.is_leaf() or (node.is_leaf() and use_chromevol_method): # Example condition
        
        node_display_name = node.name if node.name else f"Internal_Node_{i}"
        if node.is_leaf() and not node.name: node_display_name = f"Leaf_{node.get_leaves()[0].name}"


        data_row = {"node_name": node_display_name}
        
        # Common attribute 'count' (either from parsimony or ML copied over)
        data_row["inferred_chromosome_count"] = getattr(node, 'count', None)

        # ML-specific attributes
        if use_chromevol_method:
            data_row["ml_probability"] = getattr(node, 'ml_probability', None)
            # Add expected events from stochastic mapping
            for etype in StochasticMapper("","","")._summarize_event_counts({}).keys(): # Hack to get event types
                 data_row[f"expected_{etype}"] = getattr(node, f"expected_{etype}", None)
            # data_row["state_probabilities_vector"] = str(getattr(node, 'state_probabilities', [])[:10]) + "..." # Example: store part of prob vector

        if not node.is_leaf(): # Typically only save internal nodes for ancestral states report
            ancestral_states_data_list.append(data_row)
        elif use_chromevol_method and node.is_leaf(): # If ML, might want to save info for leaves too
            ancestral_states_data_list.append(data_row)


    ancestral_states_df = pd.DataFrame(ancestral_states_data_list)
    ancestral_states_df.to_csv(output_ancestral_states, index=False)
    logging.info(f"\n* Ancestral states report saved to {output_ancestral_states}")

    # --- Detailed Analysis of Key Evolutionary Transitions (using map_file) ---
    if analyze_pairs_list and map_file and os.path.exists(map_file):
        logging.info("\n--- Analysis of Key Chromosomal Rearrangement Events (from map file) ---")
        try:
            mapping_df = pd.read_csv(map_file, sep='\t')
            if len(analyze_pairs_list) % 2 != 0:
                logging.warning("--analyze_pairs requires an even number of species names (pairs). Ignoring the last unpaired species.")
                analyze_pairs_list = analyze_pairs_list[:-1]

            for i in range(0, len(analyze_pairs_list), 2):
                species_a, species_b = analyze_pairs_list[i], analyze_pairs_list[i+1]
                analyze_rearrangements(species_a, species_b, mapping_df)
        except Exception as e:
            logging.error(f"Error during synteny-based rearrangement analysis: {e}")
    elif analyze_pairs_list and not map_file:
         logging.warning("--analyze_pairs was specified, but no map_file was provided. Skipping this analysis.")


    # --- Global Quantification of Rearrangements (from tree structure) ---
    quantify_rearrangements_from_tree(tree, output_rearrangements)

    # --- Generate Annotated Tree Image ---
    plot_annotated_tree(tree, output_image, layout, show_branch_length, show_support, dpi, title="Ancestral Chromosome Evolution")


def quantify_rearrangements_from_tree(tree, output_file):
    """
    Quantifies and reports the number of fusion and fission events based on 'count' attributes
    across the entire tree. This is independent of the `map_file` synteny analysis.
    """
    logging.info("\n--- Global Analysis of Chromosome Number Evolution (from tree structure) ---")
    all_events = []
    for node in tree.traverse():
        if node.is_root() or not hasattr(node, 'count') or node.count is None \
           or not hasattr(node.up, 'count') or node.up.count is None:
            continue

        parent_count = node.up.count
        child_count = node.count
        change = child_count - parent_count

        if change == 0: continue # No change in chromosome number

        event_type = "fission" if change > 0 else "fusion"
        node_name_str = node.name if node.name else f"To_Leaf_{node.get_leaves()[0].name[:10]}"
        parent_name_str = node.up.name if node.up.name else "Root_or_Unnamed"

        all_events.append({
            "branch_to_node": node_name_str,
            "parent_node": parent_name_str,
            "event_type": event_type,
            "change_magnitude": abs(change),
            "count_from": parent_count,
            "count_to": child_count
        })

    if not all_events:
        logging.info("No chromosome number changes detected across the tree.")
        if output_file: # Create empty file or file with header
            pd.DataFrame(columns=["branch_to_node", "parent_node", "event_type", "change_magnitude", "count_from", "count_to"]).to_csv(output_file, index=False)
            logging.info(f"* Empty rearrangement events report saved to {output_file}")
        return

    events_df = pd.DataFrame(all_events)
    # Sort for consistent output, e.g., by type then magnitude
    events_df.sort_values(by=['event_type', 'change_magnitude'], ascending=[True, False], inplace=True)

    if output_file:
        events_df.to_csv(output_file, index=False)
        logging.info(f"* Tree-based rearrangement events report saved to {output_file}")

    logging.info("\nDetailed report of chromosomal changes per branch (from tree structure):")
    for _, event in events_df.iterrows():
        logging.info(f" - Branch from {event['parent_node']} to {event['branch_to_node']}: {event['event_type']} of {event['change_magnitude']} (from {event['count_from']} to {event['count_to']})")

    total_fusions = events_df[events_df['event_type'] == 'fusion']['change_magnitude'].sum()
    total_fissions = events_df[events_df['event_type'] == 'fission']['change_magnitude'].sum()

    logging.info(f"\nTotal inferred fusion events (magnitude sum): {total_fusions}")
    logging.info(f"Total inferred fission events (magnitude sum): {total_fissions}")

    logging.info("\n--- Overall Evolutionary Trend (from tree structure) ---")
    if total_fissions > total_fusions:
        logging.info("The evolutionary history of this group appears to be dominated by chromosome fissions.")
    elif total_fusions > total_fissions:
        logging.info("The evolutionary history of this group appears to be dominated by chromosome fusions.")
    else:
        logging.info("Chromosome fusions and fissions appear to be relatively balanced across the tree.")


def main():
    """Main function to parse arguments and run the analysis."""
    parser = argparse.ArgumentParser(
        description="Perform ancestral chromosome reconstruction with ChromEvol-inspired methods and visualize the results.",
        formatter_class=argparse.RawTextHelpFormatter # Allows newlines in help messages
    )

    # --- Input Files ---
    parser.add_argument('--tree', type=str, default='data/pruned_phylogeny.nwk', # Default from project structure
                        help='Path to the phylogenetic tree file in Newick format.')
    parser.add_argument('--counts', type=str, default='data/chromosome_counts.csv',
                        help='Path to the CSV file with species chromosome counts (species,chromosome_count).')
    parser.add_argument('--map', type=str, default=None, # Default to None, as it's for a specific analysis type
                        help='Path to the TSV file with chromosome mapping/synteny data (for --analyze_pairs).')

    # --- Output Files ---
    parser.add_argument('--out_image', type=str, default='results/annotated_phylogeny.svg',
                        help='Path to save the output phylogeny image. Extension determines format (svg, png, pdf).')
    parser.add_argument('--out_ancestors', type=str, default='results/ancestral_states_report.csv',
                        help='Path to save the CSV report of inferred ancestral states.')
    parser.add_argument('--out_rearrangements', type=str, default='results/rearrangement_events_report.csv',
                        help='Path to save the CSV report of rearrangement events (from tree structure).')

    # --- Analysis Parameters ---
    parser.add_argument('--root_count', type=int, default=None, # Default to None (e.g. uniform prior or estimate from data)
                        help='Constrain the ancestral chromosome number for the root of the tree (ML prior). For parsimony, this is a hard constraint.')
    parser.add_argument('--max_chr_number', type=int, default=None,
                        help='Maximum chromosome number to consider in the model. If None, inferred from data + buffer.')
    # Example for passing dict-like params for initial model values (more complex, consider config file instead)
    # parser.add_argument('--initial_model_params', type=str, help='JSON string or path to JSON file for initial model parameters.')


    # --- Plotting Parameters ---
    parser.add_argument('--layout', type=str, default='circular', choices=['circular', 'rectangular'],
                        help='Tree layout style for visualization.')
    parser.add_argument('--show_branch_length', action='store_true', default=False,
                        help='Show branch lengths in the tree visualization.')
    parser.add_argument('--show_support', action='store_true', default=False,
                        help='Show branch support values (if present in Newick) in the tree visualization.')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Dots Per Inch (DPI) for raster image formats (e.g., PNG).')

    # --- Specific Analyses ---
    parser.add_argument('--analyze_pairs', nargs='*', metavar='SPECIES', default=[],
                        help='Specify pairs of species for detailed synteny-based rearrangement analysis (requires --map file).\n'
                             'Example: --analyze_pairs SpeciesA SpeciesB SpeciesC SpeciesD')

    # --- ChromEvol Method Control ---
    parser.add_argument('--use_chromevol', action='store_true', default=False,
                        help='Use ChromEvol-inspired maximum likelihood methods. If False, uses parsimony.')
    parser.add_argument('--optimize_params', action='store_true', default=False,
                        help='Optimize model parameters using maximum likelihood (requires --use_chromevol).')
    parser.add_argument('--model_selection', action='store_true', default=False,
                        help='Perform model selection using AIC/BIC (requires --use_chromevol and implies param optimization).')

    args = parser.parse_args()

    # Create results directory if it doesn't exist
    os.makedirs("results", exist_ok=True)

    # Basic validation for file existence for required inputs at the start
    if not os.path.exists(args.tree):
        logging.error(f"Input tree file not found: {args.tree}")
        return
    if not os.path.exists(args.counts):
        logging.error(f"Input counts file not found: {args.counts}")
        return
    if args.analyze_pairs and (not args.map or not os.path.exists(args.map)):
        logging.warning(f"Map file {args.map} is required for --analyze_pairs but not found or not specified. Skipping pair analysis.")
        args.analyze_pairs = [] # Clear it so the function isn't called with missing map

    # Process initial model parameters if provided (e.g., from a JSON string or file)
    initial_model_params = None
    # if args.initial_model_params:
    #     try:
    #         if os.path.exists(args.initial_model_params):
    #             with open(args.initial_model_params, 'r') as f:
    #                 initial_model_params = json.load(f)
    #         else:
    #             initial_model_params = json.loads(args.initial_model_params)
    #     except Exception as e:
    #         logging.error(f"Error parsing --initial_model_params: {e}. Using default model parameters.")
    #         initial_model_params = None


    ancestral_reconstruction(
        tree_file=args.tree,
        counts_file=args.counts,
        map_file=args.map,
        output_image=args.out_image,
        output_ancestral_states=args.out_ancestors,
        output_rearrangements=args.out_rearrangements,
        root_constraint_count=args.root_count,
        layout=args.layout,
        show_branch_length=args.show_branch_length,
        show_support=args.show_support,
        analyze_pairs_list=args.analyze_pairs,
        dpi=args.dpi,
        use_chromevol_method=args.use_chromevol,
        optimize_model_params=args.optimize_params,
        perform_model_selection=args.model_selection,
        max_chromosome_val=args.max_chr_number, # Pass the new argument
        initial_model_params_dict=initial_model_params # Pass parsed initial params
    )


if __name__ == "__main__":
    main()
