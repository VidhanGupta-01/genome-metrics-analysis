import pytest
import pandas as pd
import sys
import os

# Make sure src is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.analysis.outlier_detection import compute_zscores, detect_outliers
from src.analysis.comparison import validate_rank, get_rank_column


# ── Sample data for testing ──────────────────────────────────

def make_sample_df():
    """Create a minimal test DataFrame mirroring Ensembl schema."""
    return pd.DataFrame({
        'assembly_id': [1, 2, 3],
        'asm_name': ['AssemblyA', 'AssemblyB', 'AssemblyC'],
        'asm_level': ['Chromosome', 'Scaffold', 'Contig'],
        'genus': ['GenusA', 'GenusB', 'GenusC'],
        'family': ['FamilyA', 'FamilyB', 'FamilyC'],
        'order_name': ['OrderA', 'OrderB', 'OrderC'],
        'contig_N50': [100000.0, 200000.0, 300000.0],
        'total_length': [1000000000.0, 2000000000.0, 3000000000.0],
    })


# ── Tests for outlier_detection ──────────────────────────────

class TestComputeZscores:

    def test_zscore_mean_is_zero(self):
        """Z-scores of a column should have mean close to 0."""
        df = make_sample_df()
        metric_cols = ['contig_N50', 'total_length']
        zscores = compute_zscores(df, metric_cols)

        for col in metric_cols:
            assert abs(zscores[col].mean()) < 1e-9, (
                f"Mean of z-scores for {col} should be ~0"
            )

    def test_zscore_constant_column_is_zero(self):
        """If all values are same, z-score should be 0 (no division by zero)."""
        df = make_sample_df()
        df['contig_N50'] = 100000.0  # constant column
        zscores = compute_zscores(df, ['contig_N50'])
        assert (zscores['contig_N50'] == 0.0).all(), (
            "Constant column should produce z-scores of 0"
        )


class TestDetectOutliers:

    def test_no_outliers_normal_data(self):
        """Normal data should produce no outliers at threshold 3.0."""
        df = make_sample_df()
        result = detect_outliers(df, threshold=3.0)
        assert result['is_outlier'].sum() == 0, (
            "Normal data should have no outliers at z > 3.0"
        )

    def test_outlier_detected_low_threshold(self):
        """With very low threshold, outliers should be detected."""
        df = make_sample_df()
        result = detect_outliers(df, threshold=0.5)
        assert result['is_outlier'].sum() > 0, (
            "Should detect outliers at very low threshold"
        )

    def test_outlier_result_has_correct_columns(self):
        """Result DataFrame should have expected columns."""
        df = make_sample_df()
        result = detect_outliers(df, threshold=3.0)
        expected_cols = {
            'assembly_id', 'asm_name', 'is_outlier',
            'flagged_metrics', 'max_zscore'
        }
        assert expected_cols.issubset(result.columns), (
            f"Missing columns in outlier result: "
            f"{expected_cols - set(result.columns)}"
        )


# ── Tests for comparison ─────────────────────────────────────

class TestValidateRank:

    def test_valid_ranks_pass(self):
        """Valid taxonomy ranks should not raise errors."""
        for rank in ['species', 'genus', 'family', 'order_name']:
            validate_rank(rank)  # should not raise

    def test_invalid_rank_raises(self):
        """Invalid rank should raise ValueError."""
        with pytest.raises(ValueError):
            validate_rank('kingdom')

    def test_rank_column_mapping(self):
        """species should map to asm_name, others to their own name."""
        assert get_rank_column('species') == 'asm_name'
        assert get_rank_column('genus') == 'genus'
        assert get_rank_column('family') == 'family'
        assert get_rank_column('order_name') == 'order_name'