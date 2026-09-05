# Historical record — not the runnable workflow

These are the original repository files, preserved verbatim for provenance.
They contain fixed lab paths, inconsistent parameters, incomplete commands,
and Markdown stored under configuration filenames. Use the root README and
`scripts/workflow.py` for the maintained entry points.

## Unresolved research provenance

- The old parameter file uses clustering 0.85 and adapter filtering 0;
  the old README describes 0.90 and 2.
- The parameter and demultiplexing files use SbfI-related settings. Library
  chemistry must be verified from the sequencing protocol before reuse.
- Clone filtering specifies a six-base oligo. Confirm actual oligo placement
  and length before applying clone filtering; it is protocol-dependent.
- The old QC script mentions removing a sample, while the newer README says
  no individuals were removed after final assembly. Neither establishes the
  complete sample history without matching logs.
- The old README reports 293,819 retained loci. The final stats file and logs
  are not present here, so that number is not presented as a verified result
  in the new overview.

Archived scripts are documentation, not validated runnable examples. The new
workflow does not change any research data or assert a corrected final analysis.
