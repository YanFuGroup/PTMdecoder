import argparse
import pandas as pd
import warnings
from pathlib import Path
from typing import Dict, List, Union

# Suppress warnings for empty slices or infinite values if they arise
warnings.filterwarnings('ignore')

def parse_metrics_from_file(filepath: Union[str, Path]) -> Dict[str, pd.Series]:
    """
    Parses the tab-separated txt file to extract stability scores: 
    'jaccard', 'support', and 'mad'.
    
    Args:
        filepath: Path to the target text file.
        
    Returns:
        A dictionary containing pandas Series for each metric.
    """
    jaccard_records = []
    support_records = []
    mad_records = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and lines starting with 'P'
            if not line or line.startswith('P'):
                continue
                
            cols = line.split('\t')
            
            # Spectrum level (Starts with 'S'): Extract jaccard
            if line.startswith('S'):
                if len(cols) >= 6:
                    for col in cols[3:6]:
                        if '=' in col:
                            key, val = col.split('=', 1)
                            if key.strip() == 'jaccard':
                                jaccard_records.append(float(val))
                                
            # Component level (Starts with other letters): Extract support and mad
            else:
                if len(cols) >= 3:
                    for col in cols[2:]:
                        if '=' in col:
                            key, val = col.split('=', 1)
                            key_clean = key.strip()
                            if key_clean == 'support':
                                support_records.append(float(val))
                            elif key_clean == 'mad':
                                mad_records.append(float(val))
                                
    return {
        'jaccard': pd.Series(jaccard_records, dtype=float),
        'support': pd.Series(support_records, dtype=float),
        'mad': pd.Series(mad_records, dtype=float)
    }

def calculate_directional_thresholds(
    metrics_dict: Dict[str, pd.Series],
    metric_quantile_map: Dict[str, List[float]]
) -> pd.DataFrame:
    """
    Calculates directional thresholds for each stability metric.

    Args:
        metrics_dict: Dictionary containing pandas Series of metrics.
        metric_quantile_map: Per-metric quantile vectors.
            Example:
            {
                'jaccard': [0.01, 0.05],
                'support': [0.01, 0.05],
                'mad': [0.99, 0.95]
            }

    Returns:
        A pandas DataFrame with metric names as index and threshold labels as columns.
    """
    results: Dict[str, Dict[str, float]] = {}

    for metric_name, series in metrics_dict.items():
        clean_series = series.dropna()
        metric_result: Dict[str, float] = {}
        quantiles = metric_quantile_map.get(metric_name, [])

        if not quantiles:
            continue

        if clean_series.empty:
            continue
        else:
            for q in quantiles:
                if q <= 0.5:
                    label = f"bottom_{int(q * 100)}%"
                else:
                    label = f"top_{int((1.0 - q) * 100)}%"
                metric_result[label] = float(clean_series.quantile(q))

        results[metric_name] = metric_result

    return pd.DataFrame.from_dict(results, orient='index')

def _validate_quantile_vector(arg_name: str, quantiles: List[float]) -> None:
    """Validates that each quantile is strictly between 0 and 1."""
    for q in quantiles:
        if not (0.0 < q < 1.0):
            raise ValueError(
                f"Invalid value in {arg_name}: {q}. Each quantile must satisfy 0 < q < 1."
            )


def main(filepath: Union[str, Path], metric_quantile_map: Dict[str, List[float]]):
    """
    Main execution pipeline for analyzing stability scores.
    """
    print(f"Reading data from: {filepath}")
    
    try:
        # 1. Parse data
        metrics = parse_metrics_from_file(filepath)
        
        print("\n--- Data Extraction Summary ---")
        for metric, series in metrics.items():
            print(f"Extracted {len(series)} values for '{metric}'")
            
        print("\n--- Active Quantile Configuration ---")
        for metric_name, quantiles in metric_quantile_map.items():
            print(f"{metric_name}: {quantiles}")

        # 2. Calculate directional thresholds with per-metric quantile vectors
        quantile_df = calculate_directional_thresholds(metrics, metric_quantile_map)
        
        # 3. Output results
        print("\n--- Stability Score Directional Thresholds ---")
        print(quantile_df.to_string())
        
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found. Please check the path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description=(
            "Analyze stability scores (Jaccard, Support, MAD) from an MS/MS report file "
            "and calculate directional thresholds: bottom 1%/5% for Jaccard and Support, "
            "top 1%/5% for MAD."
        )
    )
    
    # Add input file argument
    parser.add_argument(
        "-i", "--input", 
        type=str, 
        required=True, 
        help="Absolute or relative path to the input text file (e.g., /path/to/report_msms.txt)"
    )

    parser.add_argument(
        "--jaccard-q",
        type=float,
        nargs='+',
        default=[0.01, 0.05],
        help="Quantile vector for jaccard (e.g., --jaccard-q 0.01 0.05)"
    )

    parser.add_argument(
        "--support-q",
        type=float,
        nargs='+',
        default=[0.01, 0.05],
        help="Quantile vector for support (e.g., --support-q 0.01 0.05)"
    )

    parser.add_argument(
        "--mad-q",
        type=float,
        nargs='+',
        default=[0.99, 0.95],
        help="Quantile vector for mad (e.g., --mad-q 0.99 0.95)"
    )
    
    args = parser.parse_args()

    try:
        _validate_quantile_vector('--jaccard-q', args.jaccard_q)
        _validate_quantile_vector('--support-q', args.support_q)
        _validate_quantile_vector('--mad-q', args.mad_q)
    except ValueError as exc:
        parser.error(str(exc))

    metric_quantile_map = {
        'jaccard': args.jaccard_q,
        'support': args.support_q,
        'mad': args.mad_q
    }
    
    # Execute main function with the provided path
    main(args.input, metric_quantile_map)