#!/usr/bin/env python3
"""
Comprehensive Data Summary Generator for Chromosome Evolution Analysis
Generates publication-ready tables and statistics from the analysis results.
"""

import pandas as pd
import numpy as np
# from pathlib import Path # Not strictly needed if using os.path or direct strings
import argparse
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


def load_data(ancestral_file, events_file, counts_file):
    """Loads the required CSV files."""
    try:
        ancestral_df = pd.read_csv(ancestral_file)
        logging.info(f"Successfully loaded ancestral states from: {ancestral_file}")
    except FileNotFoundError:
        logging.error(f"Ancestral states file not found: {ancestral_file}")
        return None, None, None
    except pd.errors.EmptyDataError:
        logging.error(f"Ancestral states file is empty: {ancestral_file}")
        return None, None, None

    try:
        events_df = pd.read_csv(events_file)
        logging.info(f"Successfully loaded rearrangement events from: {events_file}")
    except FileNotFoundError:
        logging.warning(f"Rearrangement events file not found: {events_file}. Some summaries will be skipped.")
        events_df = pd.DataFrame() # Return empty DataFrame if not found
    except pd.errors.EmptyDataError:
        logging.warning(f"Rearrangement events file is empty: {events_file}. Some summaries will be skipped.")
        events_df = pd.DataFrame()


    try:
        counts_df = pd.read_csv(counts_file)
        logging.info(f"Successfully loaded chromosome counts from: {counts_file}")
    except FileNotFoundError:
        logging.error(f"Chromosome counts file not found: {counts_file}")
        return None, None, None # Counts file is essential for some stats
    except pd.errors.EmptyDataError:
        logging.error(f"Chromosome counts file is empty: {counts_file}")
        return None, None, None

    return ancestral_df, events_df, counts_df


def summarize_extant_diversity(counts_df):
    """Summarizes chromosome diversity in extant species."""
    if counts_df is None or counts_df.empty or 'chromosome_count' not in counts_df.columns:
        logging.warning("Counts data is missing or invalid. Skipping extant diversity summary.")
        return {}

    logging.info("\n1. 现存物种染色体多样性")
    logging.info("-" * 50)
    
    diversity_stats = counts_df['chromosome_count'].describe()
    num_species = len(counts_df)
    count_min = diversity_stats.get('min', np.nan)
    count_max = diversity_stats.get('max', np.nan)
    mean_count = diversity_stats.get('mean', np.nan)
    std_count = diversity_stats.get('std', np.nan)
    median_count = diversity_stats.get('50%', np.nan)

    logging.info(f"分析物种数量: {num_species}")
    logging.info(f"染色体数目范围: {count_min:.0f} - {count_max:.0f}")
    logging.info(f"平均染色体数目: {mean_count:.1f} ± {std_count:.1f}")
    logging.info(f"中位数染色体数目: {median_count:.0f}")

    summary_metrics = {
        '分析物种总数': num_species,
        '染色体数目范围': f"{count_min:.0f}-{count_max:.0f}" if pd.notna(count_min) else "N/A",
        '平均染色体数目': f"{mean_count:.1f}"  if pd.notna(mean_count) else "N/A",
    }

    if num_species >= 5:
        logging.info("\n染色体数目最多的物种:")
        top_species = counts_df.nlargest(5, 'chromosome_count')
        for _, row in top_species.iterrows():
            logging.info(f"  {row['species']}: {row['chromosome_count']} 条染色体")

        logging.info("\n染色体数目最少的物种:")
        bottom_species = counts_df.nsmallest(5, 'chromosome_count')
        for _, row in bottom_species.iterrows():
            logging.info(f"  {row['species']}: {row['chromosome_count']} 条染色体")
    return summary_metrics

def summarize_ancestral_states(ancestral_df):
    """Summarizes ancestral state reconstruction results."""
    if ancestral_df is None or ancestral_df.empty or 'node_name' not in ancestral_df.columns or 'inferred_chromosome_count' not in ancestral_df.columns:
        logging.warning("Ancestral states data is missing or invalid. Skipping ancestral states summary.")
        return
        
    logging.info("\n\n2. 祖先状态重建结果")
    logging.info("-" * 50)

    # Filter for "major" ancestral nodes (heuristic: not containing '_')
    # This might need adjustment based on actual node naming conventions from the main script.
    # A more robust way would be to have a column indicating major nodes or pass a list of such nodes.
    major_nodes = ancestral_df[~ancestral_df['node_name'].astype(str).str.contains('_', na=True)] # Handle potential NaN in node_name
    if major_nodes.empty:
        logging.info("未找到符合筛选条件的主要祖先节点 (名称不含 '_')。显示所有内部节点 (最多10个):")
        major_nodes = ancestral_df.head(10) # Show some if no "major" ones found by heuristic

    logging.info("关键祖先节点的染色体数目:")
    for _, row in major_nodes.iterrows():
        node_name = row['node_name']
        inferred_count = row['inferred_chromosome_count']
        ml_prob = row.get('ml_probability') # Use .get() for safer access

        log_msg = f"  {node_name}: {inferred_count} 条染色体"
        if pd.notna(ml_prob):
            log_msg += f" (ML概率: {float(ml_prob):.3f})"
        logging.info(log_msg)

def summarize_evolutionary_events(events_df):
    """Summarizes chromosome evolution events (fusions/fissions from tree structure)."""
    if events_df is None or events_df.empty or 'event_type' not in events_df.columns or 'change_magnitude' not in events_df.columns:
        logging.warning("Events data is missing or invalid. Skipping evolutionary events summary.")
        return {}

    logging.info("\n\n3. 染色体进化事件分析 (基于树结构)")
    logging.info("-" * 50)

    fusion_events = events_df[events_df['event_type'] == 'fusion']
    fission_events = events_df[events_df['event_type'] == 'fission']

    total_events = len(events_df)
    num_fusions = len(fusion_events)
    num_fissions = len(fission_events)

    logging.info(f"总染色体重排事件数: {total_events}")
    if total_events > 0:
        logging.info(f"  融合事件数: {num_fusions} ({num_fusions/total_events*100:.1f}%)")
        logging.info(f"  分裂事件数: {num_fissions} ({num_fissions/total_events*100:.1f}%)")
        logging.info(f"  融合:分裂 比率: {num_fusions/max(num_fissions,1):.1f}:1")

    summary_metrics = {
        '总重排事件数': total_events,
        '融合事件数': num_fusions,
        '分裂事件数': num_fissions,
        '融合:分裂 比率': f"{num_fusions/max(num_fissions,1):.1f}:1" if total_events > 0 else "N/A",
    }

    if num_fusions > 0 and 'change_magnitude' in fusion_events.columns and 'branch_to_node' in fusion_events.columns and 'count_from' in fusion_events.columns and 'count_to' in fusion_events.columns:
        major_fusions = fusion_events.nlargest(5, 'change_magnitude')
        logging.info("\n最大的融合事件:")
        for _, row in major_fusions.iterrows():
            logging.info(f"  分支至 {row['branch_to_node']}: {row['change_magnitude']} 条染色体融合 ({row['count_from']}→{row['count_to']})")
    
    if num_fissions > 0 and 'change_magnitude' in fission_events.columns and 'branch_to_node' in fission_events.columns and 'count_from' in fission_events.columns and 'count_to' in fission_events.columns:
        major_fissions = fission_events.nlargest(5, 'change_magnitude')
        logging.info("\n最大的分裂事件:")
        for _, row in major_fissions.iterrows():
            logging.info(f"  分支至 {row['branch_to_node']}: {row['change_magnitude']} 条染色体分裂 ({row['count_from']}→{row['count_to']})")
    return summary_metrics

def summarize_stochastic_mapping_rates(ancestral_df):
    """Summarizes evolutionary rates from stochastic mapping if available."""
    if ancestral_df is None or ancestral_df.empty:
        logging.warning("Ancestral states data missing. Skipping stochastic mapping summary.")
        return

    logging.info("\n\n4. 进化速率分析 (基于随机映射)")
    logging.info("-" * 50)
    
    # Dynamically find all 'expected_...' columns
    expected_event_cols = [col for col in ancestral_df.columns if col.startswith('expected_')]
    
    if not expected_event_cols:
        logging.info("在祖先状态报告中未找到随机映射的期望事件数 (expected_...) 列。")
        return

    logging.info("随机映射结果 (每个分支的期望事件数):")
    total_expected_events_summary = {}
    for col in expected_event_cols:
        event_type = col.replace('expected_', '')
        total_sum = ancestral_df[col].fillna(0).sum()
        logging.info(f"  总期望 {event_type} 事件数: {total_sum:.2f}")
        total_expected_events_summary[f'总期望 {event_type} 事件数'] = f"{total_sum:.2f}"
        
    # Calculate total activity for each branch
    ancestral_df['total_expected_activity'] = ancestral_df[expected_event_cols].fillna(0).sum(axis=1)
    
    active_branches = ancestral_df[ancestral_df['total_expected_activity'] > 0].nlargest(5, 'total_expected_activity')
    if not active_branches.empty:
        logging.info("\n进化最活跃的谱系 (基于总期望事件数):")
        for _, row in active_branches.iterrows():
            if 'node_name' in row and pd.notna(row['node_name']):
                 logging.info(f"  {row['node_name']}: {row['total_expected_activity']:.2f} 总期望事件数")
    else:
        logging.info("未发现显著的进化活跃谱系 (所有分支的期望事件总和为0或数据不足)。")
    return total_expected_events_summary


def summarize_phylogenetic_patterns(events_df):
    """Placeholder for phylogenetic pattern summary. This part is highly dataset-specific."""
    if events_df is None or events_df.empty:
        logging.warning("Events data missing. Skipping phylogenetic pattern summary.")
        return

    logging.info("\n\n5. 系统发育分布模式 (示例性)")
    logging.info("-" * 50)
    logging.info("注意: 此部分的进化枝分类是示例性的，需要根据具体研究进行定制。")
    
    # Example: Group analysis by keywords in branch names. This is a very rough heuristic.
    # A better approach would involve mapping nodes to predefined clades.
    clade_patterns = {}
    if 'branch_to_node' not in events_df.columns or 'event_type' not in events_df.columns:
        logging.warning("必要的列 (branch_to_node, event_type) 在事件数据中缺失，无法进行系统发育模式分析。")
        return

    # Simplified example: count events based on simple keywords if they appear in node names
    # This is highly dependent on how nodes are named in `events_df` (which comes from `quantify_rearrangements_from_tree`)
    keywords_example = ['Echinodermata', 'Ambulacraria', 'Deuterostomia'] # Example keywords

    for _, row in events_df.iterrows():
        branch_name = str(row['branch_to_node']) # Name of the node to which the branch leads
        event_t = row['event_type']

        assigned_clade = '其他分支' # Default
        for keyword in keywords_example:
            if keyword.lower() in branch_name.lower():
                assigned_clade = keyword
                break # Assign to first matching keyword

        clade_patterns.setdefault(assigned_clade, {'fusion': 0, 'fission': 0, 'other':0})
        if event_t in clade_patterns[assigned_clade]:
            clade_patterns[assigned_clade][event_t] += 1
        else: # Should not happen if event_type is just fusion/fission from tree structure
            clade_patterns[assigned_clade]['other'] +=1
            
    if clade_patterns:
        logging.info("主要进化枝的事件数 (基于分支名称关键词的粗略估计):")
        for clade, events_counts in clade_patterns.items():
            total_clade_events = sum(events_counts.values())
            if total_clade_events > 0:
                details = ", ".join([f"{etype}: {count}" for etype, count in events_counts.items() if count > 0])
                logging.info(f"  {clade}: {details} (总计: {total_clade_events})")
    else:
        logging.info("未能根据关键词提取进化枝模式。")


def summarize_statistics(events_df):
    """Calculates and prints various statistical summaries from event data."""
    if events_df is None or events_df.empty or 'change_magnitude' not in events_df.columns or \
       'count_from' not in events_df.columns or 'count_to' not in events_df.columns:
        logging.warning("Events data is missing or invalid for statistical summary.")
        return {}

    logging.info("\n\n6. 统计摘要")
    logging.info("-" * 50)

    changes = events_df['change_magnitude'].astype(float)
    mean_change = changes.mean()
    std_change = changes.std()
    median_change = changes.median()
    min_change = changes.min()
    max_change = changes.max()

    logging.info("重排事件的染色体数目变化幅度统计:")
    logging.info(f"  平均变化幅度: {mean_change:.1f} ± {std_change:.1f} 条染色体")
    logging.info(f"  中位数变化幅度: {median_change:.1f} 条染色体")
    logging.info(f"  变化幅度范围: {min_change:.0f} - {max_change:.0f} 条染色体")

    summary_metrics = {
        '平均变化幅度': f"{mean_change:.1f} ± {std_change:.1f}",
        '中位数变化幅度': f"{median_change:.1f}",
        '最大单次变化幅度': f"{max_change:.0f}",
    }

    initial_counts = events_df['count_from'].astype(float)
    final_counts = events_df['count_to'].astype(float)

    logging.info("\n各进化步骤的染色体数目演化:")
    logging.info(f"  变化前染色体数目范围: {initial_counts.min():.0f} - {initial_counts.max():.0f}")
    logging.info(f"  变化后染色体数目范围: {final_counts.min():.0f} - {final_counts.max():.0f}")
    return summary_metrics


def export_summary_table(summary_data_dict, output_table_file):
    """Exports the collected summary metrics to a CSV table."""
    if not summary_data_dict:
        logging.warning("无摘要数据可导出。")
        return

    logging.info("\n\n7. 生成摘要表格")
    logging.info("-" * 50)

    # Convert dictionary to DataFrame for easy CSV export
    # Ensure consistent order if Python version < 3.7 (dict order not guaranteed)
    # For newer Python, dict order is insertion order.

    # We create a list of dicts then convert to DataFrame to control order better
    df_list = [{'指标': k, '数值': v} for k, v in summary_data_dict.items()]
    summary_df = pd.DataFrame(df_list)

    try:
        summary_df.to_csv(output_table_file, index=False, encoding='utf-8-sig') # utf-8-sig for Excel compatibility with CJK
        logging.info(f"摘要表格已保存至: {output_table_file}")
    except Exception as e:
        logging.error(f"无法保存摘要表格至 {output_table_file}: {e}")


def generate_comprehensive_summary(ancestral_file, events_file, counts_file, output_table_file):
    """Main function to generate comprehensive summary tables and statistics."""

    logging.info("=" * 80)
    logging.info("染色体进化分析综合摘要")
    logging.info("=" * 80)

    ancestral_df, events_df, counts_df = load_data(ancestral_file, events_file, counts_file)

    if counts_df is None : # Counts_df is critical for some basic stats
        logging.error("由于无法加载染色体计数文件，部分摘要无法生成。")
        # return # Exit if counts_df is None as it's fundamental

    all_summary_metrics = {}

    # 1. Extant Species Diversity
    extant_metrics = summarize_extant_diversity(counts_df)
    all_summary_metrics.update(extant_metrics)

    # 2. Ancestral State Reconstruction
    summarize_ancestral_states(ancestral_df) # Mostly prints, doesn't return much for table yet

    # 3. Evolutionary Events (from tree structure)
    event_metrics = summarize_evolutionary_events(events_df)
    all_summary_metrics.update(event_metrics)

    # 4. Evolutionary Rate Analysis (Stochastic Mapping)
    stoch_map_metrics = summarize_stochastic_mapping_rates(ancestral_df) # ancestral_df should have expected_cols
    if stoch_map_metrics: all_summary_metrics.update(stoch_map_metrics)

    # 5. Phylogenetic Patterns (Example, may need customization)
    summarize_phylogenetic_patterns(events_df)

    # 6. Statistical Summary (from events_df)
    stat_metrics = summarize_statistics(events_df)
    all_summary_metrics.update(stat_metrics)

    # 7. Export Summary Table
    export_summary_table(all_summary_metrics, output_table_file)

    logging.info("\n" + "=" * 80)
    logging.info("分析摘要生成完毕！")
    logging.info("=" * 80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="生成染色体进化分析的综合摘要报告和表格。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--ancestral_states',
        type=str,
        default='results/ancestral_states_report.csv',
        help="祖先状态报告文件路径 (CSV格式)。\n默认: results/ancestral_states_report.csv"
    )
    parser.add_argument(
        '--rearrangement_events',
        type=str,
        default='results/rearrangement_events_report.csv',
        help="重排事件报告文件路径 (CSV格式)。\n默认: results/rearrangement_events_report.csv"
    )
    parser.add_argument(
        '--chromosome_counts',
        type=str,
        default='data/chromosome_counts.csv',
        help="物种染色体数目文件路径 (CSV格式)。\n默认: data/chromosome_counts.csv"
    )
    parser.add_argument(
        '--output_summary_table',
        type=str,
        default='results/analysis_summary_table.csv',
        help="输出的摘要表格文件路径 (CSV格式)。\n默认: results/analysis_summary_table.csv"
    )

    args = parser.parse_args()

    # Ensure results directory exists for output table
    os.makedirs(os.path.dirname(args.output_summary_table) or '.', exist_ok=True)

    generate_comprehensive_summary(
        ancestral_file=args.ancestral_states,
        events_file=args.rearrangement_events,
        counts_file=args.chromosome_counts,
        output_table_file=args.output_summary_table
    )
