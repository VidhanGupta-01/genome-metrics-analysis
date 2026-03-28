import pandas as pd
import matplotlib.pyplot as plt
import yaml
import os

# In production, taxonomy rank info would be fetched from NCBI Taxonomy API
# using the lowest_taxon_id field. For this prototype, taxonomy columns
# (genus, family, order_name) are pre-loaded in assembly_info.csv.

METADATA_COLS = [
    'assembly_id', 'gca_chain', 'gca_version', 'asm_type',
    'asm_level', 'asm_name', 'lowest_taxon_id', 'is_current',
    'refseq_accession', 'release_date', 'submitter',
    'genus', 'family', 'order_name'
]


def get_metric_columns(df):
    """Return only numeric metric columns, excluding metadata."""
    return [col for col in df.columns
            if col not in METADATA_COLS
            and pd.api.types.is_numeric_dtype(df[col])]


def validate_rank(rank):
    """
    Validate that the requested taxonomy rank is supported.
    Supported ranks: species (asm_name), genus, family, order_name
    """
    supported = ['species', 'genus', 'family', 'order_name']
    if rank not in supported:
        raise ValueError(
            f"Unsupported rank '{rank}'. Choose from: {supported}"
        )


def get_rank_column(rank):
    """
    Map taxonomy rank name to actual DataFrame column name.
    'species' maps to 'asm_name' since that's how Ensembl stores it.
    """
    mapping = {
        'species': 'asm_name',
        'genus': 'genus',
        'family': 'family',
        'order_name': 'order_name'
    }
    return mapping[rank]


def group_by_taxonomy(df, rank='family'):
    """
    Group assemblies by taxonomy rank.
    Returns a dict: {group_name -> subset DataFrame}

    Handles missing taxonomy gracefully by placing
    unknowns in 'Unknown' group.
    """
    validate_rank(rank)
    rank_col = get_rank_column(rank)

    # Handle missing taxonomy data
    df = df.copy()
    df[rank_col] = df[rank_col].fillna('Unknown')

    groups = {}
    for name, subset in df.groupby(rank_col):
        groups[name] = subset.reset_index(drop=True)

    print(f"\nGrouped {len(df)} assemblies into "
          f"{len(groups)} {rank}-level groups:")
    for name, subset in groups.items():
        members = ', '.join(subset['asm_name'].tolist())
        print(f"  {name}: {members}")

    return groups


def compare_metrics(df, rank='family'):
    """
    Compare metrics across taxonomy groups.
    Returns summary statistics per group per metric.
    """
    validate_rank(rank)
    rank_col = get_rank_column(rank)
    metric_cols = get_metric_columns(df)

    df = df.copy()
    df[rank_col] = df[rank_col].fillna('Unknown')

    summary = df.groupby(rank_col)[metric_cols].agg(
        ['mean', 'std', 'min', 'max', 'count']
    ).round(2)

    return summary


def print_comparison_report(summary, rank):
    """Print a readable comparison report."""
    print("\n" + "="*55)
    print(f"   TAXONOMY COMPARISON REPORT (rank: {rank})")
    print("="*55)
    print(summary.to_string())
    print("="*55)


def plot_comparison(df, rank='family', metric='contig_N50',
                    output_path=None):
    """
    Plot metric comparison across taxonomy groups as bar chart.
    Each group shows mean value with individual points overlaid.
    """
    validate_rank(rank)
    rank_col = get_rank_column(rank)

    df = df.copy()
    df[rank_col] = df[rank_col].fillna('Unknown')

    groups = df[rank_col].unique()
    colors = plt.cm.Set2.colors

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, group in enumerate(groups):
        subset = df[df[rank_col] == group]
        values = subset[metric].dropna()
        mean_val = values.mean()

        # Bar for mean value
        ax.bar(i, mean_val, color=colors[i % len(colors)],
               alpha=0.7, label=group, edgecolor='black',
               linewidth=0.5)

        # Overlay individual points
        ax.scatter(
            [i] * len(values), values,
            color='black', s=60, zorder=5
        )

        # Label each point with species name
        for j, (_, row) in enumerate(subset.iterrows()):
            if pd.notna(row[metric]):
                ax.annotate(
                    row['asm_name'],
                    (i, row[metric]),
                    textcoords="offset points",
                    xytext=(8, 4),
                    fontsize=8
                )

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=20, ha='right')
    ax.set_ylabel(metric, fontsize=11)
    ax.set_title(
        f"{metric} comparison by {rank}",
        fontsize=13, fontweight='bold'
    )
    ax.legend(title=rank.capitalize())
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150)
        print(f"Plot saved to {output_path}")

    plt.show()


if __name__ == "__main__":
    from src.db.load_data import load_and_prepare_data

    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)

    rank = config['analysis']['default_taxonomy_rank']

    df = load_and_prepare_data()

    # Group and compare
    groups = group_by_taxonomy(df, rank=rank)
    summary = compare_metrics(df, rank=rank)
    print_comparison_report(summary, rank)

    # Plot contig_N50 comparison across groups
    plot_comparison(
        df, rank=rank,
        metric='contig_N50',
        output_path="data/outputs/comparison_plot.png"
    )