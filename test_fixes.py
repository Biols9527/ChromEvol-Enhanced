#!/usr/bin/env python3

"""
Test script to validate the fixes made to ancestral_reconstruction.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
from ete3 import Tree
from ancestral_reconstruction import ChromEvolutionModel, ChromEvolLikelihoodCalculator

def test_improved_probability_normalization():
    """Test that probability normalization handles edge cases correctly"""
    print("Testing improved probability normalization...")
    
    model = ChromEvolutionModel(max_chromosome_number=10)
    
    # Test with normal branch length
    P = model.get_transition_probabilities(0.1)
    assert P.shape == (11, 11), f"Expected shape (11, 11), got {P.shape}"
    
    # Check that rows sum to 1
    row_sums = P.sum(axis=1)
    assert np.allclose(row_sums, 1.0), f"Row sums not close to 1.0: {row_sums}"
    
    # Test with zero branch length
    P_zero = model.get_transition_probabilities(0.0)
    assert np.allclose(P_zero, np.eye(11)), "Zero branch length should return identity matrix"
    
    # Test with negative branch length (should be corrected)
    P_neg = model.get_transition_probabilities(-0.1)
    assert np.allclose(P_neg, np.eye(11)), "Negative branch length should be corrected to identity matrix"
    
    print("✓ Probability normalization tests passed")

def test_missing_data_handling():
    """Test that missing data is handled consistently"""
    print("Testing missing data handling...")
    
    model = ChromEvolutionModel(max_chromosome_number=5)
    
    # Create a simple tree
    tree = Tree("(A:0.1,B:0.1);", format=1)
    
    # Data with one missing value
    data = {"A": 5, "B": None}
    
    calculator = ChromEvolLikelihoodCalculator(tree, data, model)
    
    # Test likelihood calculation with missing data
    for leaf in tree:
        if leaf.name == "A":
            likelihood = calculator.calculate_likelihood_at_node(leaf)
            # Should have likelihood 1.0 only at state 5
            expected = np.zeros(6)
            expected[5] = 1.0
            assert np.allclose(likelihood, expected), f"Expected likelihood for A: {expected}, got {likelihood}"
        elif leaf.name == "B":
            likelihood = calculator.calculate_likelihood_at_node(leaf)
            # Should have uniform likelihood for all states
            assert np.all(likelihood == 1.0), f"Expected uniform likelihood for B, got {likelihood}"
    
    print("✓ Missing data handling tests passed")

def test_parameter_validation():
    """Test that parameter validation works correctly"""
    print("Testing parameter validation...")
    
    # Test with negative parameters
    model = ChromEvolutionModel(max_chromosome_number=5)
    
    # Set negative parameters - should be corrected
    test_params = {'gain': -0.1, 'loss': 0.1, 'dupl': 0.05, 'linear_rate': 0.01}
    Q = model.build_transition_rate_matrix(test_params)
    
    # Matrix should be valid (no negative rates except diagonal)
    assert np.all(Q[Q != np.diag(Q)] >= 0), "Off-diagonal elements should be non-negative"
    
    # Test diagonal elements are negative (sum of outgoing rates)
    for i in range(model.max_chr + 1):
        assert Q[i, i] == -np.sum(Q[i, :i]) - np.sum(Q[i, i+1:]), "Diagonal elements should be negative sum of outgoing rates"
    
    print("✓ Parameter validation tests passed")

def test_thread_safety():
    """Test that likelihood calculator can be reset for thread safety"""
    print("Testing thread safety improvements...")
    
    model = ChromEvolutionModel(max_chromosome_number=5)
    tree = Tree("(A:0.1,B:0.1);", format=1)
    data = {"A": 3, "B": 4}
    
    calculator = ChromEvolLikelihoodCalculator(tree, data, model)
    
    # Calculate likelihoods
    initial_likelihood = calculator.calculate_total_log_likelihood()
    assert len(calculator.likelihoods) > 0, "Likelihoods should be cached"
    
    # Reset and recalculate
    calculator.reset_likelihoods()
    assert len(calculator.likelihoods) == 0, "Likelihoods should be cleared after reset"
    
    second_likelihood = calculator.calculate_total_log_likelihood()
    assert np.isclose(initial_likelihood, second_likelihood), "Likelihoods should be consistent after reset"
    
    print("✓ Thread safety tests passed")

def test_numerical_stability():
    """Test numerical stability improvements"""
    print("Testing numerical stability...")
    
    model = ChromEvolutionModel(max_chromosome_number=5)
    
    # Test with very large branch length
    P = model.get_transition_probabilities(1000.0)
    assert not np.any(np.isnan(P)), "Should not have NaN values even with large branch length"
    assert not np.any(np.isinf(P)), "Should not have infinite values even with large branch length"
    
    # Test with very small branch length
    P_small = model.get_transition_probabilities(1e-10)
    assert not np.any(np.isnan(P_small)), "Should not have NaN values with small branch length"
    
    print("✓ Numerical stability tests passed")

def main():
    """Run all tests"""
    print("Running tests for ancestral_reconstruction.py fixes...")
    print("=" * 60)
    
    try:
        test_improved_probability_normalization()
        test_missing_data_handling()
        test_parameter_validation()
        test_thread_safety()
        test_numerical_stability()
        
        print("=" * 60)
        print("✓ All tests passed! The fixes appear to be working correctly.")
        return 0
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())