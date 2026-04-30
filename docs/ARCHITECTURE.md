# PTMdecoder Architecture Overview

This note is for developers who need to regain context quickly. It focuses on
the current workflow entry points and class responsibilities, not release notes.

## Entry Points

- `main.m` validates command-line parameter file paths and calls `CPTMdecoder`.
- `CPTMdecoder` is the public MATLAB workflow entry point.
- `CPTMdecoderWorkflowConfig` parses a `.param` file into ordered workflow
  stages.
- `CPTMdecoderWorkflowRunner` executes enabled stages in order through a stage
  executor registry.

Most new workflow behavior should be added as a stage service and registered in
`CPTMdecoderWorkflowRunner`, rather than being wired directly into `main.m`.

## Workflow Stages

The currently supported stage names are defined in `CPTMdecoderWorkflowConfig`:

| Stage | Primary service/process | Main output role |
| --- | --- | --- |
| `msms_level` | `CMSMSLevelService` | MS/MS-level IMP discrimination and `report_msms.txt` |
| `peptide_quant` | `CPeptideQuantService` | peptide-level IMP quantification |
| `peptide_requant` | `CPeptideRequantService` | peptide-level re-quantification from checked RT ranges |
| `peptide_align_requant` | `CPeptideAlignRequantService` | alignment, aligned RT ranges, and re-quantification |
| `norm_peptide_quant` | `CNormalizationQuantService` | normalization-peptide quantification |
| `norm_peptide_requant` | `CNormalizationRequantService` | normalization-peptide re-quantification |
| `site_level` | `CSiteLevelSummary` | site-level IMP-list summary, plus uninterested records kept in peptide-level format |
| `site_level_dataset` | `CSiteLevelDatasetSummary` | site-by-dataset matrix only |
| `merge_to_pair_level` | `CMergeEachPair` | pair-level comparison |
| `merge_pairs_level` | `CMergePairs` | merged pair comparison output |

## Main Data Flow

1. `CWorkflowParamParser` reads raw parameter text into a map.
2. `CPTMdecoderWorkflowConfig.fromParamFile` converts the map into typed stage
   configs.
3. `CPTMdecoderWorkflowRunner.run` dispatches each enabled stage.
4. Stage services construct their own IO dependencies and processing executors.
5. Result objects are written through dedicated IO classes such as
   `CMS2ResultIO`, `CIMPQuantResultIO`, and `CFdrFilteredResultIO`.

## Layering

| Layer | Classes |
| --- | --- |
| Entry and workflow | `main`, `CPTMdecoder`, `CPTMdecoderWorkflowConfig`, `CPTMdecoderWorkflowRunner` |
| Stage services | `CMSMSLevelService`, `CPeptide*Service`, `CNormalization*Service` |
| Dataset IO | `CDatasetIO`, `CMgfDatasetIO`, `CMS12DatasetIO`, `CMsFileMapper`, `CFastaReader` |
| MS/MS core | `CMS2SpectrumPipeline`, `CMS2MassCalculator`, `CMS2PeakMatcher`, `CMS2QuantSolver`, `CMS2Result`, `CMS2ResultIO` |
| Peptide-level IMP quantification by XIC integration | `CIMPProcessingExecutor`, `CIMPQuantCore`, `CXIC*Utils`, `CIMPGroupAggregator`, `CIMPQuantPlotter` |
| Alignment | `CPeptideAlignRequantService`, `CIMPXICAlignRequantExecutor`, `CXICAligner`, `CXICAlign*`, `CRunAlignStrategy` implementations |
| Report models | `CIMPQuantReport`, `CIMPQuantBlock`, `CIMPQuantRecord`, `CIMPGroup`, `CIMPRawIdentManager`, `CIMPRawIdentStore` |
| Utilities | `CLogger`, `CLoggerCore`, `CPathResolver`, `CParamMapUtils`, `CStructOptionUtils`, `CModificationRegistry` |

## Current Preferred Paths

- For full runs, call `CPTMdecoder(paramFile)`.
- For adding a workflow step, add a config object or config parser, a stage
  service, and a runner registry entry.
- For MS/MS single-spectrum logic, use the `CMS2*` classes. Do not introduce new
  logic outside the current MS/MS core.
- For peptide-level quantification and re-quantification, route through
  `CIMPProcessingExecutor` from a stage service.
- For output formats, keep parsing and serialization in `*IO` classes.
