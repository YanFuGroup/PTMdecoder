"""
MaxQuant to PTMdecoder Input Converter
======================================
This script converts MaxQuant's standard `msms.txt` output into two files required 
by the PTMdecoder workflow:
1. A grouped peptide-spectrum list file (`pep_spec_file`).
2. An identification result table with detailed metrics and protein assignments.

Usage:
    python convert_mq_to_ptmdecoder.py -i <path_to_msms.txt> -o <path_to_pep_spec.txt> -r <path_to_result.txt>
"""

import pandas as pd
import argparse
import sys
from typing import Tuple

# ------------------------------------------------------------------------
# Module 1: Core Parsing & Calculation Logic
# ------------------------------------------------------------------------

def parse_modified_sequence(mod_seq: str) -> Tuple[str, str]:
    """
    Parses MaxQuant's 'Modified sequence' to extract modifications and their locations.
    Rule: Left '_' is pos 0. Each capital letter adds 1 to pos. Right '_' adds 1 to pos (L+1).
    Handles nested parentheses.
    """
    if pd.isna(mod_seq) or mod_seq == '' or mod_seq == '_':
        return '-', '-'
    
    mods, locs = [], []
    current_pos, i = 0, 0
    
    while i < len(mod_seq):
        char = mod_seq[i]
        
        if char == '_':
            if i == 0: 
                current_pos = 0  # N-terminal underscore, position 0
            else:
                current_pos += 1 # C-terminal underscore, position becomes total amino acid count + 1
            i += 1
            
        elif char == '(':
            j = i + 1
            open_brackets = 1
            while j < len(mod_seq) and open_brackets > 0:
                if mod_seq[j] == '(': open_brackets += 1
                elif mod_seq[j] == ')': open_brackets -= 1
                j += 1
            
            mod_name = mod_seq[i+1:j-1]
            mods.append(mod_name)
            locs.append(str(current_pos))
            i = j
            
        else:
            current_pos += 1     # Normal amino acid, position + 1
            i += 1

    if not mods:
        return '-', '-'
    return ','.join(mods), ','.join(locs)


def calculate_precursor_mass(df: pd.DataFrame) -> pd.Series:
    """Calculates precursor neutral mass from m/z and Charge."""
    try:
        # Use errors='raise' to strictly enforce numeric conversion. 
        # If dirty data (e.g., letters, corrupted text) is encountered, a ValueError is raised immediately.
        mz_arr = pd.to_numeric(df['m/z'], errors='raise')
        charge_arr = pd.to_numeric(df['Charge'], errors='raise')
        
        mass_proton = 1.007276
        precursor_mass = (mz_arr - mass_proton) * charge_arr
        
        # Since all values are guaranteed to be valid numbers at this point, 
        # we can directly format them to 6 decimal places without checking pd.notna().
        return precursor_mass.map(lambda x: f"{x:.6f}")
        
    except ValueError as e:
        # Specifically catch errors related to failed numeric conversions.
        print(f"Error: Dirty data detected during mass calculation! Unconvertible values found in 'm/z' or 'Charge' columns.")
        print(f"Detailed error message: {e}")
        sys.exit(1)  # Exit the Python script immediately with status code 1.
        
    except Exception as e:
        # Catch any other unexpected errors.
        print(f"Error: A critical error occurred while calculating precursor_neutral_mass: {e}")
        sys.exit(1)

# ------------------------------------------------------------------------
# Module 2: Data I/O and Validation
# ------------------------------------------------------------------------

def load_and_validate_data(input_file: str, strict_con_filter: bool = False) -> pd.DataFrame:
    """Loads MaxQuant msms.txt, validates required columns, and cleans data."""
    print(f"Reading MaxQuant results from: {input_file} ...")
    try:
        df = pd.read_csv(input_file, sep='\t', dtype=str)
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found. Please check the path.")
        sys.exit(1)

    required_cols = [
        'Raw file', 'Scan number', 'Charge', 'Sequence', 
        'Mass', 'Mass error [Da]', 'm/z', 'Number of matches', 
        'Proteins', 'Modified sequence', 'Score', 'Reverse'
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in input file: {missing_cols}")
        sys.exit(1)
        
    # 0. Filter out Reverse hits (Decoy database)
    # MaxQuant uses '+' for reverse hits. We only keep rows where 'Reverse' is NaN or empty string.
    initial_count = len(df)
    df = df[df['Reverse'].isna() | (df['Reverse'].str.strip() == '')]
    filtered_count = len(df)
    
    if initial_count - filtered_count > 0:
        print(f"[Info] Filtered out {initial_count - filtered_count} reverse (decoy) hits.")
        
    # 1. Drop strictly invalid rows for basic identification metrics
    df = df.dropna(subset=['Raw file', 'Scan number', 'Charge', 'Sequence'])
    
    # 2. Strictly validate Proteins column: Raise error if NA is found!
    if df['Proteins'].isna().any():
        # Extract all rows where 'Proteins' is missing
        missing_proteins_df = df[df['Proteins'].isna()]
        
        print("Error: Missing protein assignment information detected! The 'Proteins' column contains null values (NA), which is unacceptable in the current workflow.")
        print("Please check your MaxQuant search settings or the FASTA database file.")
        print(f"A total of {len(missing_proteins_df)} rows are missing 'Proteins'.")
        print("Below are the first few locations of the erroneous data (line numbers correspond to the original txt file):")
        
        # Iterate through the first 5 erroneous rows to print details without flooding the console
        for idx, row in missing_proteins_df.head(5).iterrows():
            # Pandas index starts at 0, and the header takes up the 1st line, so physical line number = idx + 2
            line_num = idx + 2 
            raw_name = row['Raw file']
            scan_num = row['Scan number']
            print(f"  -> Line {line_num} | Raw file: {raw_name} | Scan number: {scan_num}")
            
        sys.exit(1)

    # 3. Fill missing values in 'Mass error [Da]' and 'Number of matches' with '0' to ensure consistent output format.
    df['Mass error [Da]'] = df['Mass error [Da]'].fillna('0')
    df['Number of matches'] = df['Number of matches'].fillna('0')
    
    # 4. Format Proteins column for downstream compatibility
    # Replace semicolons with commas for multiple protein assignments.
    # Filtering strategy for proteins starting with 'CON__':
    # - lenient mode (default): remove only CON__ entries in-place; drop rows only if nothing remains.
    # - strict mode: drop the entire row if any CON__ entry is present.
    df['Proteins'] = df['Proteins'].astype(str).str.replace(';', ',', regex=False)

    def _contains_con_protein(protein_str: str) -> bool:
        if protein_str is None:
            return False
        parts = [p.strip() for p in protein_str.split(',') if p.strip() != '']
        return any(p.startswith('CON__') for p in parts)

    def _filter_con_proteins(protein_str: str) -> str:
        if protein_str is None:
            return ''
        parts = [p.strip() for p in protein_str.split(',') if p.strip() != '']
        filtered = [p for p in parts if not p.startswith('CON__')]
        return ','.join(filtered)

    if strict_con_filter:
        strict_drop_mask = df['Proteins'].apply(_contains_con_protein)
        if strict_drop_mask.any():
            num_removed = strict_drop_mask.sum()
            print(f"[Info] Strict mode enabled: removing {num_removed} rows containing at least one 'CON__' protein assignment.")
            df = df[~strict_drop_mask]
    else:
        df['Proteins'] = df['Proteins'].apply(_filter_con_proteins)

        # Drop rows where Proteins is empty after filtering out CON__ entries
        empty_mask = df['Proteins'].astype(str).str.strip() == ''
        if empty_mask.any():
            num_removed = empty_mask.sum()
            print(f"[Info] Removing {num_removed} rows where only 'CON__' proteins were assigned (now empty after filtering).")
            df = df[~empty_mask]
    
    return df

# ------------------------------------------------------------------------
# Module 3: Output Generators
# ------------------------------------------------------------------------

def generate_pep_spec_list(df: pd.DataFrame, output_path: str, dataset_suffix: str = ''):
    """Generates the grouped peptide-spectrum list file (Task 1)."""
    print("Generating grouped peptide-spectrum list...")
    
    # 1. Filter out rows where modification is '-'
    # Call the existing parse_modified_sequence function to parse 'Modified sequence'.
    # Keep the row if the parsed modification part is not '-'.
    is_modified = df['Modified sequence'].apply(lambda x: parse_modified_sequence(x)[0] != '-')
    filtered_df = df[is_modified]

    # 2. Group the filtered DataFrame by Sequence and write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        grouped = filtered_df.groupby('Sequence', sort=False)
        for sequence, group in grouped:
            f.write(f"{sequence}\n")
            for _, row in group.iterrows():
                raw = row['Raw file']
                scan = row['Scan number']
                charge = row['Charge']
                col1 = f"{raw}{dataset_suffix}.mgf"
                col2 = scan
                f.write(f"{col1}\t{col2}\n")
    print(f"[Success] Peptide-spectrum list saved to: {output_path}")


def generate_ident_result_table(df: pd.DataFrame, output_path: str, dataset_suffix: str = ''):
    """Parses modifications and generates the detailed identification result table (Task 2)."""
    print("Parsing modifications and computing masses for identification result table...")
    
    precursor_mass_str = calculate_precursor_mass(df)
    parsed_mods = df['Modified sequence'].apply(parse_modified_sequence)
    modifications = [x[0] for x in parsed_mods]
    mod_locations = [x[1] for x in parsed_mods]

    df_result = pd.DataFrame({
        'Site': '-', 
        'DatasetName': df['Raw file'].astype(str) + dataset_suffix + '.mgf',
        'Scan': df['Scan number'],
        'Spectrum': df['Scan number'],
        'Charge': df['Charge'],
        'Calc_neutral_pepmass': df['Mass'],
        'precursor_neutral_mass': precursor_mass_str,
        'massdiff': df['Mass error [Da]'],
        'num_match_ions': df['Number of matches'],
        'peptide': df['Sequence'],
        'protein': df['Proteins'],
        'modification': modifications,
        'modificationlocation': mod_locations,
        'Score': df['Score']
    })

    df_result.to_csv(output_path, sep='\t', index=False)
    print(f"[Success] Identification result table saved to: {output_path}")

# ------------------------------------------------------------------------
# Main Orchestrator
# ------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert MaxQuant msms.txt to PTMdecoder input formats."
    )
    parser.add_argument('-i', '--input', required=True, help="Path to input MaxQuant msms.txt")
    parser.add_argument('-o', '--output', required=True, help="Path to output pepSpecFile.txt")
    parser.add_argument('-r', '--result', required=True, help="Path to output filtered_result.txt")
    parser.add_argument('-s', '--suffix', default='', help="Optional dataset suffix to insert between Raw file and '.mgf' (supports empty string)")
    parser.add_argument(
        '-c', '--strict-con-filter',
        action='store_true',
        help="Enable strict contaminant filtering: drop any row containing a protein that starts with 'CON__'. Default is lenient mode.",
    )
    
    args = parser.parse_args()
    
    # 1. Load and Validate
    df = load_and_validate_data(args.input, strict_con_filter=args.strict_con_filter)
    
    # 2. Generate Task 1 File
    generate_pep_spec_list(df, args.output, args.suffix)
    
    # 3. Generate Task 2 File
    generate_ident_result_table(df, args.result, args.suffix)
    
    print("All tasks completed! You are ready to run PTMdecoder.")

if __name__ == "__main__":
    main()