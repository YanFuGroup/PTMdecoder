import pytest
import pandas as pd
import numpy as np

# Import functions from the new modularized analyzer script
from stability_score_analyzer import parse_metrics_from_file, calculate_directional_thresholds

@pytest.fixture
def mock_txt_file(tmp_path):
    """
    Creates a temporary text file with mock data mimicking the target format.
    """
    file_content = """P	meta1	meta2	meta3
S	1	2	jaccard=0.90	vif_all=1.2	vif_reported=1.1
S	1	2	jaccard=0.10	vif_all=2.5	vif_reported=2.0
S	1	2	jaccard=0.50	vif_all=1.0	vif_reported=1.0
M	A	B	support=100	mad=0.05	vif=1.5
X	C	D	support=10	mad=0.95	vif=2.0
Y	E	F	support=50	mad=0.50	vif=1.1
"""
    file_path = tmp_path / "mock_report_msms.txt"
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(file_content)
        
    return file_path

def test_parse_metrics_from_file(mock_txt_file):
    """
    Tests if the parsing function correctly extracts the stability metrics.
    """
    metrics = parse_metrics_from_file(mock_txt_file)
    
    # Check if all keys exist
    assert "jaccard" in metrics
    assert "support" in metrics
    assert "mad" in metrics
    
    # Check the lengths of the extracted series
    assert len(metrics["jaccard"]) == 3
    assert len(metrics["support"]) == 3
    assert len(metrics["mad"]) == 3
    
    # Check specific extracted values
    assert list(metrics["jaccard"]) == [0.90, 0.10, 0.50]
    assert list(metrics["support"]) == [100.0, 10.0, 50.0]
    assert list(metrics["mad"]) == [0.05, 0.95, 0.50]

def test_calculate_directional_thresholds():
    """
    Tests if the quantile calculation logic works correctly on known data.
    """
    # Create mock series (1 to 100 makes percentiles easy to verify)
    mock_data = {
        'jaccard': pd.Series(np.arange(1, 101, dtype=float)),  # 1.0 to 100.0
        'support': pd.Series(np.arange(1, 101, dtype=float)), 
        'mad': pd.Series(np.arange(1, 101, dtype=float))
    }
    
    metric_quantile_map = {
        'jaccard': [0.01, 0.05],
        'support': [0.01, 0.05],
        'mad': [0.99, 0.95]
    }
    result_df = calculate_directional_thresholds(mock_data, metric_quantile_map)

    # Check index (metrics)
    assert "jaccard" in result_df.index
    assert "support" in result_df.index
    assert "mad" in result_df.index

    # Check directional label columns
    assert "bottom_1%" in result_df.columns
    assert "bottom_5%" in result_df.columns
    assert "top_1%" in result_df.columns
    assert "top_5%" in result_df.columns

    # In an array of 1 to 100, the 0.01 quantile is approximately 1.99
    assert np.isclose(result_df.loc['jaccard', 'bottom_1%'], 1.99)
    assert np.isclose(result_df.loc['jaccard', 'bottom_5%'], 5.95)
    assert np.isclose(result_df.loc['mad', 'top_1%'], 99.01)
    assert np.isclose(result_df.loc['mad', 'top_5%'], 95.05)

def test_empty_file_handling(tmp_path):
    """
    Tests if the script handles completely empty files without crashing.
    """
    empty_file = tmp_path / "empty.txt"
    empty_file.write_text("")
    
    metrics = parse_metrics_from_file(empty_file)
    metric_quantile_map = {
        'jaccard': [0.01, 0.05],
        'support': [0.01, 0.05],
        'mad': [0.99, 0.95]
    }
    result_df = calculate_directional_thresholds(metrics, metric_quantile_map)
    
    # Empty input should not crash and should produce an empty result table.
    assert len(metrics["jaccard"]) == 0
    assert result_df.empty