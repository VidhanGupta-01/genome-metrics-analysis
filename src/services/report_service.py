import pandas as pd
import yaml
from src.db.load_data import load_and_prepare_data
from src.analysis.outlier_detection import detect_outliers
from src.analysis.comparison import group_by_taxonomy, get_metric_columns


def load_config():
    with open('config.yaml', 'r') as f:
        return yaml.safe_load(f)


def get_assembly_stats(df, assembly_id):
    """
    Get all metrics for a specific assembly.
    Returns a Series with metric name → value.
    """
    row = df[df['assembly_id'] == assembly_id]

    if row.empty:
        raise ValueError(f"Assembly ID {assembly_id} not found in data.")

    metric_cols = get_metric_columns(df)
    return row[metric_cols].iloc[0]


def get_metric_percentiles(df, assembly_id):
    """
    For each metric, show where this assembly ranks
    compared to all others (as percentile).

    Example: contig_N50 at 80th percentile means
    this genome has better contig_N50 than 80% of others.
    """
    metric_cols = get_metric_columns(df)
    row = df[df['assembly_id'] == assembly_id][metric_cols].iloc[0]

    percentiles = {}
    for col in metric_cols:
        all_values = df[col].dropna()
        value = row[col]
        pct = (all_values < value).sum() / len(all_values) * 100
        percentiles[col] = round(pct, 1)

    return percentiles


def get_taxonomy_context(df, assembly_id, rank):
    """
    Show how this assembly compares within its taxonomy group.
    Returns group name and group summary stats.
    """
    from src.analysis.comparison import get_rank_column
    rank_col = get_rank_column(rank)

    row = df[df['assembly_id'] == assembly_id].iloc[0]
    group_name = row[rank_col]

    # Get all assemblies in same group
    group_df = df[df[rank_col] == group_name]
    metric_cols = get_metric_columns(df)

    summary = group_df[metric_cols].agg(['mean', 'min', 'max']).round(2)
    return group_name, group_df['asm_name'].tolist(), summary


def generate_report(assembly_id):
    """
    Main function: generate a full report for a given assembly.
    """
    config = load_config()
    threshold = config['analysis']['outlier_threshold']
    rank = config['analysis']['default_taxonomy_rank']

    df = load_and_prepare_data()

    # Basic info
    row = df[df['assembly_id'] == assembly_id].iloc[0]

    print("\n" + "="*60)
    print("         PER-GENOME ANNOTATION METRICS REPORT")
    print("="*60)
    print(f"  Assembly Name   : {row['asm_name']}")
    print(f"  Assembly ID     : {assembly_id}")
    print(f"  Assembly Type   : {row['asm_type']}")
    print(f"  Assembly Level  : {row['asm_level']}")
    print(f"  Submitter       : {row['submitter']}")
    print(f"  Release Date    : {row['release_date']}")
    print(f"  Is Current      : {row['is_current']}")
    print("="*60)

    # Metric values
    print("\n  METRICS:")
    stats = get_assembly_stats(df, assembly_id)
    for metric, value in stats.items():
        print(f"    {metric:<25} {value:>15,.0f}")

    # Percentile ranks
    print("\n  PERCENTILE RANKS (vs all assemblies):")
    percentiles = get_metric_percentiles(df, assembly_id)
    for metric, pct in percentiles.items():
        bar = "█" * int(pct / 10) + "░" * (10 - int(pct / 10))
        print(f"    {metric:<25} {bar}  {pct:.0f}th percentile")

    # Outlier check
    print("\n  OUTLIER CHECK:")
    outlier_df = detect_outliers(df, threshold=threshold)
    row_outlier = outlier_df[outlier_df['assembly_id'] == assembly_id].iloc[0]

    if row_outlier['is_outlier']:
        print(f"    ⚠ FLAGGED as outlier (threshold: z > {threshold})")
        print(f"    Flagged metrics: {row_outlier['flagged_metrics']}")
    else:
        print(f"    ✓ No outliers detected (threshold: z > {threshold})")
        print(f"    Max z-score: {row_outlier['max_zscore']}")

    # Taxonomy context
    print(f"\n  TAXONOMY CONTEXT (rank: {rank}):")
    group_name, members, summary = get_taxonomy_context(
        df, assembly_id, rank
    )
    print(f"    Group     : {group_name}")
    print(f"    Members   : {', '.join(members)}")
    print(f"\n    Group stats:")
    print(summary.to_string())

    print("\n" + "="*60)


if __name__ == "__main__":
    # Generate report for Human genome (assembly_id = 1)
    generate_report(assembly_id=1)