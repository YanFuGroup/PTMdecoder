# PTMdecoder Maintainer Guide

This guide keeps the class layout readable as the codebase grows.

## Documentation Rules

- Keep `README.md` focused on user installation and running instructions.
- Put developer-facing explanations in `docs/`.
- When adding a public workflow stage, update both `docs/ARCHITECTURE.md` and
  `docs/CLASS_GUIDE.md`.
- When changing a file format, update the relevant `*IO` class comments and any
  docs that mention the output file.

## Class Ownership Rules

- Workflow dispatch belongs in `CPTMdecoderWorkflowRunner`.
- Parameter-to-stage construction belongs in `CPTMdecoderWorkflowConfig` and
  stage config classes.
- File parsing and serialization belong in `*IO` classes.
- Dataset access belongs in dataset IO classes, not in algorithms.
- Stage services may compose dependencies, but algorithm details should stay in
  core classes, executors, or utility classes.
- Logging should go through `CLogger`.
- Path joining and output directory creation should go through `CPathResolver`
  when possible.

## Adding a Workflow Stage

1. Add or reuse a config class that validates the stage settings.
2. Add a stage name constant in `CPTMdecoderWorkflowConfig`.
3. Build the stage from `fromParamFile`.
4. Implement a `*Service` or process class with a `run()` method.
5. Register the stage in `CPTMdecoderWorkflowRunner.buildStageExecutorRegistry`.
6. Add focused tests for config parsing and stage dispatch.
7. Update `docs/ARCHITECTURE.md` and `docs/CLASS_GUIDE.md`.

## Adding an Algorithm Helper

Before adding a new helper class, check whether the function belongs to an
existing owner:

- MS/MS mass, peak matching, and solver logic: `CMS2*`.
- XIC preprocessing, peak selection, and area logic: `CXIC*Utils`.
- Peptide-level quant/re-quant execution: `CIMPProcessingExecutor` or
  `CIMPQuantCore`.
- Alignment sub-steps: `CXICAlign*`.

Prefer small static utility methods for stateless math and signal operations.
Prefer data model methods only when they operate on model invariants.

## Test Expectations

- Config classes should have tests for defaults, required fields, and invalid
  values.
- IO classes should have round-trip tests when practical.
- Solver and signal-processing changes should include numerical edge cases.
- Stage services should be tested with mocks or small fixtures rather than full
  large datasets.

## Review Checklist

- Is the new class name consistent with an existing family?
- Is there one clear owner for each new responsibility?
- Did the change avoid duplicating parsing, path, or logging logic?
- Are output file names and formats documented near their IO classes?
- Are the developer docs updated when the workflow surface changes?
