import pandas as pd
import yaml
import os

def load_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(os.path.dirname(__file__), '../../config.yaml')
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_raw_metrics(raw_csv_path):
    """
    Load key-value format metrics CSV.
    Expected columns: assembly_id, metrics_name, metrics_value
    """
    df = pd.read_csv(raw_csv_path)

    required_cols = {'assembly_id', 'metrics_name', 'metrics_value'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"raw.csv must contain columns: {required_cols}")

    return df

def pivot_metrics(df):
    """
    Convert key-value format to wide format for analysis.
    
    Input:
        assembly_id | metrics_name | metrics_value
        1           | contig_N50   | 125000
    
    Output:
        assembly_id | contig_N50 | total_length | ...
        1           | 125000     | 2800000000   | ...
    """
    wide_df = df.pivot_table(
        index='assembly_id',
        columns='metrics_name',
        values='metrics_value',
        aggfunc='first'
    ).reset_index()

    wide_df.columns.name = None
    return wide_df

def load_assembly_info(assembly_csv_path):
    """Load assembly metadata (gca_chain, asm_type, asm_level etc.)"""
    return pd.read_csv(assembly_csv_path)

def load_and_prepare_data():
    """
    Main function: loads, pivots, and merges all data.
    Returns a clean wide-format DataFrame ready for analysis.
    """
    config = load_config()

    raw_path = config['data']['raw_metrics']
    assembly_path = config['data']['assembly_info']

    raw_df = load_raw_metrics(raw_path)
    wide_df = pivot_metrics(raw_df)
    assembly_df = load_assembly_info(assembly_path)

    wide_df['assembly_id'] = wide_df['assembly_id'].astype(int)
    assembly_df['assembly_id'] = assembly_df['assembly_id'].astype(int)
    merged_df = pd.merge(wide_df, assembly_df, on='assembly_id', how='left')

    print(f"Loaded {len(merged_df)} assemblies with {len(merged_df.columns)} fields")
    return merged_df

if __name__ == "__main__":
    df = load_and_prepare_data()
    print(df.head())