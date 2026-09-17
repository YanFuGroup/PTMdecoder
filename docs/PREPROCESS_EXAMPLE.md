# Mascot preprocess example for analysis callers

PTMdecoder provides reusable building blocks for converting Mascot result data
into the inputs consumed by its MS/MS workflow. It intentionally does not own a
dataset-specific preprocess pipeline. Paths, run names, publication rules, and
scientific selectors belong to the calling analysis repository.

This document is the stable, documentation-only example. A runnable minimal
example may be added later, but it must use 3-5 fictional PSMs and must not
turn this skeleton back into a dataset-specific template.

## Ownership boundary

PTMdecoder owns:

- Mascot DAT reading with `ReadDatResult` and `ReadDatResultFolder`.
- Score sorting and target/decoy/group classification with `JudgeGroup`.
- The existing GF/SF/TF calculations in `ComputeFDR`.
- The 14-column filtered-result table format through
  `CFdrFilteredResultIO` and `write_mascot_result_table`.
- The peptide/spectrum text format through
  `write_peptide_spectra_list_file`.

The calling analysis repository owns:

- Input and output paths.
- Run names and batch orchestration.
- `TagType`, decoy/group tags, the FDR threshold, and which GF/SF/TF result is
  used.
- All pre-FDR chemical or experimental filters.
- PSM uniqueness policy.
- All post-FDR selectors, including peptide, charge, mass, matched-ion, and
  modification rules.
- Final sorting, output publication, overwrite protection, and archival.

Do not add project-specific paths, run lists, selectors, or publication rules
back to PTMdecoder. Keep the customized orchestration in the analysis repository.

## Minimal caller skeleton

The following skeleton is intentionally not runnable until every placeholder is
replaced by the calling project. Copy the completed orchestration into that
project; do not submit the customized copy to PTMdecoder.

```matlab
ptmdecoderRoot = '<path-to-PTMdecoder>';
preprocessDir = fullfile(ptmdecoderRoot, 'FDR_control_generate_pep_spec_list');
addpath(preprocessDir);

mascotDatRoot = '<folder-containing-run-subdirectories-with-.dat-files>';
outputRoot = '<new-explicit-output-directory>';
if isfolder(outputRoot)
    error('Refusing to overwrite the existing output directory: %s', outputRoot);
end
mkdir(outputRoot);

% Every value below is an analysis-policy decision.
tagType = '<Protein-or-Modification>';
decoyTag = '<your-decoy-tag>';
groupTag = {'<your-group-tag>'};
fdrThreshold = <your-fdr-threshold>;
selectedFdrMode = '<GF-or-SF-or-TF>';

validFdrModes = {'GF', 'SF', 'TF'};
if ~any(strcmp(selectedFdrMode, validFdrModes))
    error('selectedFdrMode must be GF, SF, or TF.');
end

% PTMdecoder reads all .dat files in each run folder of mascotDatRoot.
result = ReadDatResultFolder(mascotDatRoot);

% Apply analysis-owned pre-FDR rules here, before JudgeGroup and ComputeFDR.
% Leave this step out entirely when the analysis has no such rule.

[DecoyType, GroupType, ~, scores, numrst, I] = JudgeGroup( ...
    result, tagType, decoyTag, groupTag);
[FDR, Iid, threshold, finalFDR] = ComputeFDR( ...
    DecoyType, GroupType, scores, numrst, I, fdrThreshold);

% The caller explicitly chooses one of the computed modes. PTMdecoder does not
% provide a project default and does not restrict callers to SF.
groupedResult = result(I(~GroupType));
filteredResult = result(Iid.(selectedFdrMode));

write_mascot_result_table( ...
    groupedResult, fullfile(outputRoot, 'group_result_mascot.txt'));
write_mascot_result_table( ...
    filteredResult, fullfile(outputRoot, 'filtered_result_mascot.txt'));

% Define and test PSM uniqueness and all post-FDR selectors in the calling
% project. Historical analysis workflows removed every PSM sharing an ambiguous
% (DatasetName, Scan) key, but PTMdecoder does not prescribe that policy here.
selectedResult = filteredResult;  % Replace with caller-owned selections.

% The pepSpec writer preserves caller order and expects equal peptides to be
% contiguous. Sort in the caller before writing.
[~, sortedIndex] = sort({selectedResult.peptide});
write_peptide_spectra_list_file( ...
    selectedResult(sortedIndex), fullfile(outputRoot, 'pepSpecFile.txt'));
```

The `FDR` output contains all three modes. Inspect `FDR.GF`, `FDR.SF`, and
`FDR.TF` when evaluating the policy; `Iid`, `threshold`, and `finalFDR` have the
same three fields.

## Compatibility notes

- `ComputeFDR` retains the legacy GF/SF/TF behavior and calls `robustfit`, so
  this skeleton requires MATLAB and the Statistics Toolbox.
- `CFdrFilteredResultIO` reads and writes the current 14-column text table.
  Numeric values are text after a read/write round trip.
- `write_peptide_spectra_list_file` does not sort, deduplicate, or validate
  peptide grouping. The caller must make those decisions before writing.
- `ReadDatResult` preserves the current Mascot DAT parser assumptions. Projects
  with unusual DAT layouts should characterize their files in the analysis
  repository before relying on the reader.
- A runnable demo, when added, must use only fictional PSMs and write to an
  automatically cleaned temporary directory. It must not use research data or
  values derived from research data.
