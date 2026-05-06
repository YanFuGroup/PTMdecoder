# PTMdecoder Class Guide

Use this guide as a "which class do I need?" index.

## By Task

| Task | Start with | Usually also involved |
| --- | --- | --- |
| Run a complete parameter file | `CPTMdecoder` | `CPTMdecoderWorkflowConfig`, `CPTMdecoderWorkflowRunner` |
| Add or reorder workflow stages | `CPTMdecoderWorkflowConfig`, `CPTMdecoderWorkflowRunner` | one new `*Service` or process class |
| Parse workflow parameters | `CWorkflowParamParser` | `CPTMdecoderWorkflowParamKeys`, `CParamMapUtils`, stage `*Config` classes |
| MS/MS-level IMP discrimination | `CMSMSLevelService` | `CMS2SpectrumPipeline`, `CMS2MassCalculator`, `CMS2PeakMatcher`, `CMS2QuantSolver` |
| Read/write MS/MS results | `CMS2ResultIO` | `CMS2Result` |
| Peptide-level quantification | `CPeptideQuantService` | `CPeptideRawIdentAssembler`, `CIMPProcessingExecutor`, `CIMPQuantResultIO` |
| Peptide-level re-quantification | `CPeptideRequantService` | `CIMPQuantResultIO`, `CIMPQuantReport`, `CIMPProcessingExecutor` |
| Alignment plus re-quantification | `CPeptideAlignRequantService` | `CIMPXICAlignRequantExecutor`, `CXICAligner`, `CRunAlignStrategy` |
| XIC preprocessing or peak detection | `CXICPreprocessUtils`, `CXICSignalUtils` | `CXICPeakUtils`, `CXICAreaUtils` |
| XIC peak area and final ratio logic | `CXICAreaUtils`, `CXICPeakUtils` | `CIMPQuantCore` |
| IMP group aggregation | `CIMPGroupAggregator` | `CIMPGroup`, `CIMPProcessingExecutor` |
| Draw XIC plots | `CXICDrawService`, `CIMPQuantPlotter` | `CXICDrawServiceConfig` |
| Site-level summaries | `CSiteLevelSummary`, `CSiteLevelDatasetSummary` | `CSiteLevel*Config` |
| Pairwise comparison output | `CMergeEachPair`, `CMergePairs` | `CMergeEachPairConfig`, `CMergePairsConfig` |
| Logging | `CLogger` | `CLoggerCore` |
| Path handling | `CPathResolver` | stage config output path fields |
| Modification masses | `CModificationRegistry` | `CConstant` |

## Main Class Families

### Workflow and Stage Services

- `CPTMdecoder` owns top-level execution.
- `CPTMdecoderWorkflowConfig` owns stage construction from parameter files.
- `CPTMdecoderWorkflowRunner` owns dispatch and stage error boundaries.
- `CMSMSLevelService`, `CPeptideQuantService`, `CPeptideRequantService`,
  `CPeptideAlignRequantService`, `CNormalizationQuantService`, and
  `CNormalizationRequantService` are stage services.

Stage services should be thin orchestration classes. They may create IO objects,
configs, and executors, but detailed algorithms should live in core classes or
utilities.

### Data Access

- `CMgfDatasetIO` reads MGF spectra by dataset and spectrum name, and builds
  cached TITLE/scan-number to precursor-charge maps.
- `CMS12DatasetIO` reads MS1/MS2 files and supports XIC-related access.
- `CMsFileMapper` maps MGF names to MS1/MS2 stems.
- `CFastaReader` reads FASTA sequences.
- `CPepProtService` resolves peptide-to-protein context.

### MS/MS Core

- `CMS2SpectrumPipeline` coordinates one-spectrum processing.
- `CMS2MassCalculator` builds theoretical masses and ion m/z values.
- `CMS2PeakMatcher` matches theoretical and experimental peaks.
- `CMS2QuantSolver` solves IMP abundance and diagnostics.
- `CMS2Result` and `CMS2ResultIO` hold and serialize stage output.

### Peptide-Level IMP Quantification and XIC Integration

- `CPeptideRawIdentAssembler` converts MS/MS results or FDR entries into
  `CIMPRawIdentManager` objects.
- `CIMPProcessingExecutor` is the executor used by quant, re-quant, and drawing
  services.
- `CIMPQuantCore` contains core IMP group quantification logic.
- `CXICSignalUtils`, `CXICPreprocessUtils`, `CXICPeakUtils`, and
  `CXICAreaUtils` contain XIC-specific pure or mostly stateless helpers.
- `CIMPQuantReport`, `CIMPQuantBlock`, and `CIMPQuantRecord` are the report data
  model.
- `CIMPQuantReport` also provides report-derived lookup maps, including peptide
  RT-range maps used by re-quantification.
- `CIMPGroupAggregator` groups IMP candidates within each raw file by mass and
  charge. Each aggregated `CIMPGroup` is then processed by
  `CIMPProcessingExecutor` for quantification, re-quantification, or drawing.

### Alignment

- `CPeptideAlignRequantService` orchestrates quant-before-alignment,
  alignment, report writing, and aligned re-quantification.
- `CIMPXICAlignRequantExecutor` builds aligned RT range maps.
- `CXICAligner` aligns IMPs between runs.
- `CXICAlignAnchorSelector`, `CXICAlignTransform`, `CXICAlignScore`, and
  `CXICAlignPeakSelector` implement alignment sub-steps.
- `PairwiseRunAlignStrategy` and `ReferenceRunAlignStrategy` decide run pairs.

## Naming Hints

- `*Service`: stage-level orchestration.
- `*Config`: validated settings for a service, executor, or summary process.
- `*IO`: file format parsing and writing.
- `*Utils`: stateless helpers or low-state algorithm utilities.
- `*Executor`: reusable execution engine called by services.
- `*Report`, `*Block`, `*Record`, `*Result`, `*Store`, `*Manager`: data models.

## When Adding a New Class

Add a class only when one of these is true:

- A new workflow stage needs a clear orchestration boundary.
- A file format needs dedicated parsing/writing logic.
- Shared algorithm code is reused by multiple services or tests.
- A data structure has invariants that are worth naming and testing.

Otherwise, prefer adding a method to the existing owner class.
