import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Same metadata columns as pca.py
METADATA_COLS = [
    'assembly_id', 'gca_chain', 'gca_version', 'asm_type',
    'asm_level', 'asm_name', 'lowest_taxon_id', 'is_current',
    'refseq_accession', 'release_date', 'submitter'
]


def get_metric_columns(df):
    """Return only numeric metric columns, excluding metadata."""
    return [col for col in df.columns
            if col not in METADATA_COLS
            and pd.api.types.is_numeric_dtype(df[col])]


def compute_zscores(df, metric_cols):
    """
    Compute z-scores for each metric column.
    Z-score = (value - mean) / std

    Returns a DataFrame of z-scores with same shape as input metrics.
    """
    zscores = df[metric_cols].copy()
    for col in metric_cols:
        mean = df[col].mean()
        std = df[col].std()
        if std == 0:
            # If all values are same, z-score is 0 (no variation)
            zscores[col] = 0.0
        else:
            zscores[col] = (df[col] - mean) / std
    return zscores


def detect_outliers(df, threshold=3.0):
    """
    Detect outlier genomes based on z-score threshold.

    For each genome, checks every metric.
    If any metric has |z-score| > threshold, genome is flagged.

    Returns a DataFrame with:
    - assembly_id, asm_name
    - is_outlier (True/False)
    - flagged_metrics (which metrics triggered the flag)
    - max_zscore (worst z-score across all metrics)
    """
    metric_cols = get_metric_columns(df)
    zscores = compute_zscores(df, metric_cols)

    results = []

    for i, row in df.iterrows():
        zscore_row = zscores.loc[i]

        # Find which metrics exceed threshold
        flagged = {
            col: round(zscore_row[col], 2)
            for col in metric_cols
            if abs(zscore_row[col]) > threshold
        }

        max_z = zscore_row.abs().max()

        results.append({
            'assembly_id': row['assembly_id'],
            'asm_name': row['asm_name'],
            'asm_level': row['asm_level'],
            'is_outlier': len(flagged) > 0,
            'flagged_metrics': flagged,
            'max_zscore': round(max_z, 2)
        })

    return pd.DataFrame(results)


def print_outlier_report(outlier_df):
    """
    Print a human-readable outlier report.
    """
    print("\n" + "="*55)
    print("       GENOME OUTLIER DETECTION REPORT")
    print("="*55)

    outliers = outlier_df[outlier_df['is_outlier'] == True]

    if outliers.empty:
        print("No outliers detected with current threshold.")
    else:
        print(f"Found {len(outliers)} outlier(s):\n")
        for _, row in outliers.iterrows():
            print(f"  Assembly : {row['asm_name']}")
            print(f"  Level    : {row['asm_level']}")
            print(f"  Max Z    : {row['max_zscore']}")
            print(f"  Flagged  : {row['flagged_metrics']}")
            print()

    print("="*55)
    print("All assemblies summary:\n")
    print(outlier_df[['asm_name', 'is_outlier', 'max_zscore']].to_string(index=False))
    print("="*55)


def plot_outliers(df, outlier_df, output_path=None):
    """
    Plot z-scores as a heatmap-style bar chart per genome.
    Highlights which metrics are extreme for each assembly.
    """
    metric_cols = get_metric_columns(df)
    zscores = compute_zscores(df, metric_cols)
    zscores['asm_name'] = df['asm_name'].values

    fig, ax = plt.subplots(figsize=(11, 5))

    x = range(len(metric_cols))
    width = 0.15
    colors = plt.cm.Set2.colors

    for i, (_, row) in enumerate(zscores.iterrows()):
        values = [row[col] for col in metric_cols]
        positions = [pos + i * width for pos in x]
        ax.bar(positions, values, width=width,
               label=row['asm_name'], color=colors[i % len(colors)])

    # Draw threshold lines
    ax.axhline(3, color='red', linestyle='--',
               linewidth=1, label='Threshold (+3)')
    ax.axhline(-3, color='red', linestyle='--', linewidth=1)

    ax.set_xticks([pos + width * 2 for pos in x])
    ax.set_xticklabels(metric_cols, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel("Z-Score", fontsize=11)
    ax.set_title("Z-Score per Metric per Genome", fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=8, title="Assembly")
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150)
        print(f"Plot saved to {output_path}")

    plt.show()


if __name__ == "__main__":
    from src.db.load_data import load_and_prepare_data
    import yaml

    # Load threshold from config
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    threshold = config['analysis']['outlier_threshold']

    df = load_and_prepare_data()
    outlier_df = detect_outliers(df, threshold=threshold)
    print_outlier_report(outlier_df)
    plot_outliers(df, outlier_df,
                  output_path="data/outputs/outlier_plot.png")