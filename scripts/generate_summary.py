#!/usr/bin/env python3
"""
Comprehensive Data Summary Generator for Chromosome Evolution Analysis
Generates publication-ready tables and statistics from the analysis results.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def generate_comprehensive_summary():
    """Generate comprehensive summary tables and statistics"""
    
    print("=" * 80)
    print("COMPREHENSIVE CHROMOSOME EVOLUTION ANALYSIS SUMMARY")
    print("=" * 80)
    
    # Load the analysis results
    try:
        ancestral_df = pd.read_csv('comprehensive_ancestral_states.csv')
        events_df = pd.read_csv('comprehensive_rearrangements.csv')
        counts_df = pd.read_csv('chromosome_counts.csv')
    except FileNotFoundError:
        print("Warning: Using default result files...")
        ancestral_df = pd.read_csv('ancestral_states_report.csv')
        events_df = pd.read_csv('rearrangement_events_report.csv')
        counts_df = pd.read_csv('chromosome_counts.csv')
    
    # 1. SPECIES CHROMOSOME DIVERSITY TABLE
    print("\n1. EXTANT SPECIES CHROMOSOME DIVERSITY")
    print("-" * 50)
    
    diversity_stats = counts_df['chromosome_count'].describe()
    print(f"Number of species analyzed: {len(counts_df)}")
    print(f"Chromosome count range: {diversity_stats['min']:.0f} - {diversity_stats['max']:.0f}")
    print(f"Mean chromosome count: {diversity_stats['mean']:.1f} ± {diversity_stats['std']:.1f}")
    print(f"Median chromosome count: {diversity_stats['50%']:.0f}")
    
    # Top and bottom species by chromosome count
    print("\nSpecies with highest chromosome counts:")
    top_species = counts_df.nlargest(5, 'chromosome_count')
    for _, row in top_species.iterrows():
        print(f"  {row['species']}: {row['chromosome_count']} chromosomes")
        
    print("\nSpecies with lowest chromosome counts:")
    bottom_species = counts_df.nsmallest(5, 'chromosome_count')
    for _, row in bottom_species.iterrows():
        print(f"  {row['species']}: {row['chromosome_count']} chromosomes")
    
    # 2. ANCESTRAL STATE RECONSTRUCTION SUMMARY
    print("\n\n2. ANCESTRAL STATE RECONSTRUCTION RESULTS")
    print("-" * 50)
    
    # Filter for major ancestral nodes
    major_nodes = ancestral_df[~ancestral_df['node_name'].str.contains('_', na=False)]
    print("Key ancestral chromosome numbers:")
    for _, row in major_nodes.iterrows():
        ml_prob = row.get('ml_probability', 'N/A')
        if ml_prob != 'N/A' and pd.notna(ml_prob):
            print(f"  {row['node_name']}: {row['inferred_chromosome_count']} chromosomes (ML prob: {float(ml_prob):.3f})")
        else:
            print(f"  {row['node_name']}: {row['inferred_chromosome_count']} chromosomes")
    
    # 3. EVOLUTIONARY EVENTS ANALYSIS
    print("\n\n3. CHROMOSOME EVOLUTION EVENTS ANALYSIS")
    print("-" * 50)
    
    # Count event types
    fusion_events = events_df[events_df['type'] == 'fusion']
    fission_events = events_df[events_df['type'] == 'fission']
    
    print(f"Total chromosomal rearrangement events: {len(events_df)}")
    print(f"  Fusion events: {len(fusion_events)} ({len(fusion_events)/len(events_df)*100:.1f}%)")
    print(f"  Fission events: {len(fission_events)} ({len(fission_events)/len(events_df)*100:.1f}%)")
    print(f"  Fusion:Fission ratio: {len(fusion_events)/max(len(fission_events),1):.1f}:1")
    
    # Major fusion events
    if len(fusion_events) > 0:
        major_fusions = fusion_events.nlargest(5, 'change')
        print("\nLargest fusion events:")
        for _, row in major_fusions.iterrows():
            print(f"  {row['branch']}: {row['change']} chromosomes fused ({row['from_count']}→{row['to_count']})")
    
    # Major fission events  
    if len(fission_events) > 0:
        major_fissions = fission_events.nlargest(5, 'change')
        print("\nLargest fission events:")
        for _, row in major_fissions.iterrows():
            print(f"  {row['branch']}: {row['change']} chromosomes split ({row['from_count']}→{row['to_count']})")
    
    # 4. EVOLUTIONARY RATE ANALYSIS
    print("\n\n4. EVOLUTIONARY RATE ANALYSIS")
    print("-" * 50)
    
    # Check if ML columns are available
    ml_columns = ['expected_gains', 'expected_losses', 'expected_duplications']
    if all(col in ancestral_df.columns for col in ml_columns):
        print("Stochastic mapping results (expected events per branch):")
        
        # Calculate total expected events
        gains_total = ancestral_df['expected_gains'].fillna(0).sum()
        losses_total = ancestral_df['expected_losses'].fillna(0).sum()
        dupl_total = ancestral_df['expected_duplications'].fillna(0).sum()
        
        print(f"  Total expected gains: {gains_total:.2f}")
        print(f"  Total expected losses: {losses_total:.2f}")
        print(f"  Total expected duplications: {dupl_total:.2f}")
        
        # Most active branches
        ancestral_df['total_activity'] = (ancestral_df['expected_gains'].fillna(0) + 
                                        ancestral_df['expected_losses'].fillna(0) + 
                                        ancestral_df['expected_duplications'].fillna(0))
        
        active_branches = ancestral_df.nlargest(5, 'total_activity')
        if len(active_branches) > 0:
            print("\nMost evolutionarily active lineages:")
            for _, row in active_branches.iterrows():
                if row['total_activity'] > 0:
                    print(f"  {row['node_name']}: {row['total_activity']:.2f} expected events")
    else:
        print("Maximum likelihood rate analysis not available in this dataset.")
    
    # 5. PHYLOGENETIC PATTERNS
    print("\n\n5. PHYLOGENETIC DISTRIBUTION PATTERNS")
    print("-" * 50)
    
    # Group analysis by major clades
    clade_patterns = {}
    for _, row in events_df.iterrows():
        branch = row['branch']
        if 'Echinodermata' in branch:
            clade = 'Echinodermata'
        elif 'Ambulacraria' in branch:
            clade = 'Ambulacraria'
        elif 'Deuterostomia' in branch:
            clade = 'Deuterostomia'
        else:
            clade = 'Terminal branches'
            
        if clade not in clade_patterns:
            clade_patterns[clade] = {'fusion': 0, 'fission': 0}
        clade_patterns[clade][row['type']] += 1
    
    print("Events per major clade:")
    for clade, events in clade_patterns.items():
        total = events['fusion'] + events['fission']
        if total > 0:
            print(f"  {clade}: {events['fusion']} fusions, {events['fission']} fissions (total: {total})")
    
    # 6. STATISTICAL SUMMARY
    print("\n\n6. STATISTICAL SUMMARY")
    print("-" * 50)
    
    # Change magnitude analysis
    changes = events_df['change'].astype(float)
    print(f"Rearrangement magnitude statistics:")
    print(f"  Mean change: {changes.mean():.1f} ± {changes.std():.1f} chromosomes")
    print(f"  Median change: {changes.median():.1f} chromosomes")
    print(f"  Range: {changes.min():.0f} - {changes.max():.0f} chromosomes")
    
    # Chromosome number evolution
    initial_counts = events_df['from_count'].astype(float)
    final_counts = events_df['to_count'].astype(float)
    
    print(f"\nChromosome number evolution:")
    print(f"  Initial state range: {initial_counts.min():.0f} - {initial_counts.max():.0f}")
    print(f"  Final state range: {final_counts.min():.0f} - {final_counts.max():.0f}")
    
    # 7. EXPORT SUMMARY TABLE
    print("\n\n7. GENERATING SUMMARY TABLES")
    print("-" * 50)
    
    # Create comprehensive summary table
    summary_table = pd.DataFrame({
        'Metric': [
            'Total species analyzed',
            'Chromosome count range',
            'Mean chromosome count',
            'Total rearrangement events',
            'Fusion events',
            'Fission events', 
            'Fusion:Fission ratio',
            'Mean change magnitude',
            'Largest single change'
        ],
        'Value': [
            len(counts_df),
            f"{counts_df['chromosome_count'].min()}-{counts_df['chromosome_count'].max()}",
            f"{counts_df['chromosome_count'].mean():.1f}",
            len(events_df),
            len(fusion_events),
            len(fission_events),
            f"{len(fusion_events)/max(len(fission_events),1):.1f}:1",
            f"{changes.mean():.1f}",
            f"{changes.max():.0f}"
        ]
    })
    
    # Save summary table
    summary_table.to_csv('analysis_summary_table.csv', index=False)
    print("Summary table saved to: analysis_summary_table.csv")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE - All summary data generated successfully!")
    print("=" * 80)

if __name__ == "__main__":
    generate_comprehensive_summary()
