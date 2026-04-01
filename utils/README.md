# PTMdecoder Utilities

This directory contains utility scripts to support data preparation, format conversion, and other auxiliary tasks for the PTMdecoder workflow.

## Available Scripts

### 1. MaxQuant Format Converter
* **Script:** `convert_mq_to_ptmdecoder.py`
* **Description:** Converts MaxQuant's standard `msms.txt` output into the specific formats required by PTMdecoder. It generates two essential files:
  1. **Grouped peptide-spectrum list** (`pepSpecFile.txt`): Formatted for RT alignment and core quantification.
  2. **Identification result table** (`filtered_result_mascot.txt`): A detailed TSV file with calculated precursor masses, dynamic modification parsing, and proper PTM positioning mappings.

* **Usage:**
  Run the script from the command line, providing the input `msms.txt` and the desired output paths for the two generated files:

  ```bash
  python convert_mq_to_ptmdecoder.py \
      -i /path/to/MaxQuant/txt/msms.txt \
      -o /path/to/output/pepSpecFile.txt \
      -r /path/to/output/filtered_result_mascot.txt
  ```

* **Advanced PTM Parsing Note:**
  This script features a state-machine parser for the `Modified sequence` column from MaxQuant. It intelligently handles nested brackets (e.g., `Acetyl (Protein N-term)`) and precisely calculates modification positions based on the rule that the N-terminus `_` is at position 0, and subsequent amino acids increment the position index.

---
*Note: For detailed documentation, internal logic, and modular functions of each tool, please refer to the docstrings within the respective Python files.*