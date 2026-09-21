# Runnable Mascot preprocess example for analysis callers

PTMdecoder provides reusable building blocks for converting Mascot result data
into the inputs consumed by its MS/MS workflow. It intentionally does not own a
dataset-specific preprocess pipeline. Paths, run names, publication rules, and
scientific selectors belong to the calling analysis repository.

The repository includes a runnable teaching example at
[`run_mascot_preprocess_example.m`](../examples/mascot_preprocess/run_mascot_preprocess_example.m).
It uses only fictional data, writes to an automatically cleaned temporary
directory, and does not define scientific defaults or a stable public API.

## Run the example

From MATLAB, point `ptmdecoderRoot` at this source checkout and run:

```matlab
ptmdecoderRoot = '<path-to-PTMdecoder>';
exampleDir = fullfile(ptmdecoderRoot, 'examples', 'mascot_preprocess');
addpath(exampleDir);
exampleResult = run_mascot_preprocess_example();
```

The function locates the repository from its own file, temporarily adds the
repository root and `FDR_control_generate_pep_spec_list` to the MATLAB path,
and restores the original path before returning. It prints every generated
file before deleting its temporary output directory. The same text remains
available under `exampleResult.outputText`.

The example reads the public fixture
[`minimal_mascot.dat`](../examples/mascot_preprocess/fixtures/minimal_mascot.dat).
This is a small parser-oriented subset containing 12 fictional PSMs. It is not
a complete example of the Mascot DAT specification and must not be treated as
research data.

## What the example demonstrates

The runnable call chain is:

```text
ReadDatResult -> JudgeGroup -> ComputeFDR
              -> write_mascot_result_table
              -> write_peptide_spectra_list_file
```

The fixture includes six in-group targets, four in-group decoys, one
out-of-group target, and one out-of-group decoy. This is enough to exercise the
existing transfer-FDR fit instead of its small-sample `[0, 0]` fallback.

Every policy value in the example is deliberately fictional and exists only to
make the output easy to inspect:

```matlab
tagType = 'Protein';
decoyTag = 'DECOY_';
groupTag = 'EXAMPLE_GROUP';
fdrThreshold = 0.5;
selectedFdrMode = 'GF';
```

These values are not recommended defaults. In particular, the threshold and
selected GF mode have no scientific meaning. The example applies no pre-FDR
filter, PSM uniqueness policy, or post-FDR selector. It sorts only by peptide
before calling the pepSpec writer because that writer expects equal peptides
to be contiguous.

The temporary directory contains three files while the example is running:

| File | Example contents |
| --- | --- |
| `group_result_mascot.txt` | 14-column table with the 10 in-group target and decoy PSMs |
| `filtered_result_mascot.txt` | 14-column table with the four in-group targets passing example-only GF filtering |
| `pepSpecFile.txt` | peptide-grouped spectrum list derived from those four targets |

`exampleResult` exposes the example-only policy, counts, FDR diagnostics, the
three text snapshots, and `cleanupConfirmed`. It does not return temporary file
paths that no longer exist.

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

Do not copy the example's fictional values into an analysis without replacing
and validating them. Do not add project-specific paths, run lists, selectors,
or publication rules back to PTMdecoder. Keep customized orchestration in the
analysis repository.

## Adapting the call chain

Use the runnable function as the maintained source example, then implement the
following analysis-owned decisions in the calling repository:

1. Choose a single DAT file for `ReadDatResult`, or a non-recursive run folder
   for `ReadDatResultFolder`.
2. Apply any justified pre-FDR rules before `JudgeGroup` and `ComputeFDR`.
3. Supply and test the analysis-specific tags, threshold, and selected FDR
   mode. `ComputeFDR` returns `GF`, `SF`, and `TF` fields in `FDR`, `Iid`,
   `threshold`, and `finalFDR`; PTMdecoder does not choose a project default.
4. Define uniqueness and post-FDR selectors in the analysis repository.
5. Sort equal peptides contiguously before writing pepSpec output.
6. Own publication paths, overwrite protection, and archival outside
   PTMdecoder.

For a single protein-group tag, pass a character vector or string scalar such
as `'EXAMPLE_GROUP'`. A one-element cell is not a valid single-tag argument for
the current `JudgeGroup` implementation.

## Compatibility notes

- The known source baseline is MATLAB R2022a. `ComputeFDR` calls `robustfit`,
  so running this source example also requires Statistics and Machine Learning
  Toolbox. Installing MATLAB or MATLAB Runtime R2022a alone does not guarantee
  that this toolbox is installed and licensed.
- The example uses platform-neutral MATLAB path and temporary-directory APIs,
  but it does not create a new Linux or macOS support commitment. The compiled
  application's supported platform remains documented in the repository
  README.
- `CFdrFilteredResultIO` writes the current 14-column tab-delimited table, but
  its reader splits rows on either tabs or runs of spaces. Numeric values become
  text after a round trip. Fields containing spaces can shift columns, and
  tokens after the first 14 are silently ignored. This remains a legacy
  compatibility limitation rather than a whitespace-safe TSV contract.
- `write_peptide_spectra_list_file` does not sort, deduplicate, or validate
  peptide grouping. The caller must make those decisions before writing.
- `ReadDatResult` preserves the current Mascot parser assumptions. Projects
  with unusual DAT layouts must characterize their files in the analysis
  repository before relying on the reader.

## Non-goals

The example is not a production runner, general framework, dataset template,
or new stable API. It does not restore removed dataset-specific workflows and
does not implement batch processing, publishing, migration, promotion,
archival, or scientific selection policy.
