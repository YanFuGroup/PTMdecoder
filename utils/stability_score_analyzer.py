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

def calculate_quantiles(metrics_dict: Dict[str, pd.Series], quantiles: List[float]) -> pd.DataFrame:
    """
    Calculates specified quantiles for the given stability metrics.
    
    Args:
        metrics_dict: Dictionary containing pandas Series of metrics.
        quantiles: List of float values representing the quantiles to calculate (e.g., [0.01, 0.05]).
        
    Returns:
        A pandas DataFrame with metrics as columns and quantiles as index.
    """
    results = {}
    for metric_name, series in metrics_dict.items():
        if series.dropna().empty:
            # Handle empty data gracefully
            results[metric_name] = {q: float('nan') for q in quantiles}
        else:
            results[metric_name] = series.quantile(quantiles).to_dict()
            
    return pd.DataFrame(results)

def main(filepath: Union[str, Path]):
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
            
        # 2. Calculate quantiles
        target_quantiles = [0.01, 0.05]
        quantile_df = calculate_quantiles(metrics, target_quantiles)
        
        # 3. Output results
        print("\n--- Stability Score Quantiles ---")
        print(quantile_df.to_string())
        
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found. Please check the path.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Analyze stability scores (Jaccard, Support, MAD) from an MS/MS report file and calculate their 0.01 and 0.05 quantiles."
    )
    
    # Add input file argument
    parser.add_argument(
        "-i", "--input", 
        type=str, 
        required=True, 
        help="Absolute or relative path to the input text file (e.g., /path/to/report_msms.txt)"
    )
    
    args = parser.parse_args()
    
    # Execute main function with the provided path
    main(args.input)