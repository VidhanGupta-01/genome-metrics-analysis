import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os

# Columns that are metadata, not metrics
METADATA_COLS = [
    'assembly_id', 'gca_chain', 'gca_version', 'asm_type',
    'asm_level', 'asm_name', 'lowest_taxon_id', 'is_current',
    'refseq_accession', 'release_date', 'submitter'
]

def get_metric_columns(df):
    """
    Returns only the numeric metric columns,
    excluding all metadata columns.
    """
    return [col for col in df.columns
            if col not in METADATA_COLS
            and pd.api.types.is_numeric_dtype(df[col])]


def standardize_metrics(df, metric_cols):
    """
    Standardize metrics using StandardScaler.
    Each metric will have mean=0 and std=1.
    This ensures no single metric dominates PCA.
    """
    scaler = StandardScaler()
    scaled = scaler.fit_transform(df[metric_cols])
    return pd.DataFrame(scaled, columns=metric_cols)


def run_pca(df, n_components=2):
    """
    Run PCA on standardized metrics.
    Returns DataFrame with PC1, PC2 and assembly metadata.
    """
    metric_cols = get_metric_columns(df)

    if len(metric_cols) < 2:
        raise ValueError("Need at least 2 numeric metric columns for PCA")

    # Standardize
    scaled_df = standardize_metrics(df, metric_cols)

    # Run PCA
    pca = PCA(n_components=n_components)
    components = pca.fit_transform(scaled_df)

    # Build result DataFrame
    result_df = pd.DataFrame(
        components,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )

    # Add metadata back
    result_df['assembly_id'] = df['assembly_id'].values
    result_df['asm_name'] = df['asm_name'].values
    result_df['asm_level'] = df['asm_level'].values

    # Explained variance
    explained = pca.explained_variance_ratio_ * 100
    print(f"PC1 explains {explained[0]:.1f}% variance")
    print(f"PC2 explains {explained[1]:.1f}% variance")

    return result_df, explained


def plot_pca(result_df, explained, output_path=None):
    """
    Plot PCA scatter plot with species labels.
    """
    fig, ax = plt.subplots(figsize=(9, 6))

    # Color by asm_level
    levels = result_df['asm_level'].unique()
    colors = plt.cm.Set2.colors

    for i, level in enumerate(levels):
        subset = result_df[result_df['asm_level'] == level]
        ax.scatter(
            subset['PC1'],
            subset['PC2'],
            label=level,
            color=colors[i % len(colors)],
            s=100,
            edgecolors='black',
            linewidths=0.5
        )

    # Label each point with species name
    for _, row in result_df.iterrows():
        ax.annotate(
            row['asm_name'],
            (row['PC1'], row['PC2']),
            textcoords="offset points",
            xytext=(8, 4),
            fontsize=9
        )

    ax.set_xlabel(f"PC1 ({explained[0]:.1f}% variance)", fontsize=11)
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}% variance)", fontsize=11)
    ax.set_title("PCA of Genome Assembly Metrics", fontsize=13, fontweight='bold')
    ax.legend(title="Assembly Level")
    ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=150)
        print(f"Plot saved to {output_path}")

    plt.show()


if __name__ == "__main__":
    from src.db.load_data import load_and_prepare_data

    df = load_and_prepare_data()
    result_df, explained = run_pca(df)
    plot_pca(result_df, explained, output_path="data/outputs/pca_plot.png")