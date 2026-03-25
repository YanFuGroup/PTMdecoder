import pandas as pd
import pytest
import sys

from convert_MaxQuant_to_PTMdecoder import (
    parse_modified_sequence,
    calculate_precursor_mass,
    load_and_validate_data,
    generate_pep_spec_list,
    generate_ident_result_table
)

# ==========================================
# Module 1: Unit Tests (Test accuracy of core logic)
# ==========================================

def test_parse_modified_sequence():
    """Test if parsing of various modifications strictly follows the rules."""
    # 1. No modification
    assert parse_modified_sequence("_PEPTIDE_") == ("-", "-")
    
    # 2. N-terminal modification only (Position 0)
    assert parse_modified_sequence("_(Acetyl (Protein N-term))PEPTIDE_") == ("Acetyl (Protein N-term)", "0")
    
    # 3. Internal modification (Count capital letters)
    assert parse_modified_sequence("_PEP(Phospho (STY))TIDE_") == ("Phospho (STY)", "3")
    
    # 4. C-terminal modification only (Length is 7, C-term should be 7+1=8)
    assert parse_modified_sequence("_PEPTIDE_(Amidation)") == ("Amidation", "8")
    
    # 5. Extremely complex mixture: N-term, nested internal brackets, and C-term coexist
    assert parse_modified_sequence("_(Ac)PE(Phospho (ST))PTIDE_(Amidation)") == ("Ac,Phospho (ST),Amidation", "0,2,8")
    
    # 6. Null value or empty exception handling
    assert parse_modified_sequence(pd.NA) == ("-", "-")
    assert parse_modified_sequence("_") == ("-", "-")

def test_calculate_precursor_mass():
    """Test mass calculation formula: (m/z - 1.007276) * Charge."""
    df = pd.DataFrame({'m/z': ['400.0'], 'Charge': ['2']})
    res = calculate_precursor_mass(df)
    expected_mass = (400.0 - 1.007276) * 2
    
    # Verify if exactly 6 decimal places are kept
    assert res.iloc[0] == f"{expected_mass:.6f}"

# ==========================================
# Module 2: Exception Tests (Test error interception)
# ==========================================

def test_missing_protein_raises_exit(tmp_path):
    """Test if the program strictly intercepts and exits when 'Proteins' is missing in forward hits."""
    input_file = tmp_path / "dirty_msms.txt"
    
    # Forge dirty data missing 'Proteins'
    data = {
        'Raw file': ['RawA'], 'Scan number': ['100'], 'Charge': ['2'], 'Sequence': ['PEP'],
        'Mass': ['800'], 'Mass error [Da]': ['0'], 'm/z': ['400'], 'Number of matches': ['10'],
        'Proteins': [float('nan')], 'Modified sequence': ['_PEP_'], 'Score': ['100'], 'Reverse': ['']
    }
    pd.DataFrame(data).to_csv(input_file, sep='\t', index=False)

    # Catch SystemExit exception
    with pytest.raises(SystemExit) as excinfo:
        load_and_validate_data(str(input_file))
    
    # Verify if exit code is 1 (indicates error exit)
    assert excinfo.value.code == 1

# ==========================================
# Module 3: End-to-End Integration Tests 
# ==========================================

def test_full_pipeline_with_decoys(tmp_path):
    """Test the full pipeline of reading files, filtering decoys, and correctly generating output files."""
    # tmp_path is a clean temporary directory provided by pytest, destroyed automatically after tests
    input_file = tmp_path / "msms.txt"
    pep_out = tmp_path / "pep_spec.txt"
    res_out = tmp_path / "result.txt"

    # 1. Forge comprehensive test data including decoys (Reverse) and modifications
    data = {
        'Raw file': ['RawA', 'RawB', 'RawC'],
        'Scan number': ['100', '200', '300'],
        'Charge': ['2', '3', '2'],
        'Sequence': ['PEPTIDE', 'SEQUENCE', 'DECOY'],
        'Mass': ['800.0', '900.0', '1000.0'],
        'Mass error [Da]': ['0.01', '-0.02', '0'],
        'm/z': ['401.0', '301.0', '501.0'],
        'Number of matches': ['10', '15', '5'],
        'Proteins': ['ProtA', 'ProtB', 'REV_ProtC'],
        'Modified sequence': ['_PEPTIDE_', '_(Ac (N-term))SEQUEN(Phos (S))CE_', '_DECOY_'],
        'Score': ['100', '120', '20'],
        'Reverse': ['', float('nan'), '+'] # 1st row empty string, 2nd row NaN, 3rd row is decoy (+)
    }
    pd.DataFrame(data).to_csv(input_file, sep='\t', index=False)

    # 2. Run data loading (Test filtering mechanism)
    df_valid = load_and_validate_data(str(input_file))
    
    # Assert: Originally 3 rows, should be exactly 2 rows left after filtering decoys
    assert len(df_valid) == 2
    assert 'REV_ProtC' not in df_valid['Proteins'].values

    # 3. Run file generation modules
    generate_pep_spec_list(df_valid, str(pep_out))
    generate_ident_result_table(df_valid, str(res_out))

    # 4. Assert if pep_spec_list file meets MATLAB requirements
    with open(pep_out, 'r', encoding='utf-8') as f:
        lines = f.read().splitlines()
        
        # Verify grouping headers
        assert "PEPTIDE" in lines
        assert "SEQUENCE" in lines
        
        # Verify concatenated format of mgf and dta
        assert "RawA.mgf\tRawA.100.100.2.0.dta" in lines

    # 5. Assert if Result table content is parsed correctly
    df_res = pd.read_csv(res_out, sep='\t', dtype=str)
    assert len(df_res) == 2
    
    # Check the first result (no modification)
    assert df_res.iloc[0]['Site'] == '-'
    assert df_res.iloc[0]['modification'] == '-'
    
    # Check the second result (with modification) for correct mapping
    assert df_res.iloc[1]['modification'] == "Ac (N-term),Phos (S)"
    assert df_res.iloc[1]['modificationlocation'] == "0,6"
    assert df_res.iloc[1]['Spectrum'] == "RawB.200.200.3.0.dta"