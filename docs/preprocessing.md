# Before assembly

The maintained workflow begins with **demultiplexed, single-end reads**, one
sample per input file. It does not automatically demultiplex raw sequencing pools.

1. Confirm the library's restriction enzyme, barcode position, and whether a
   random oligo was incorporated for identifying PCR clones.
2. Use Stacks `clone_filter` only when appropriate for that library design.
   Determine whether filtering should precede or follow demultiplexing from
   whether oligos are unique to a sample or an entire library.
3. Configure `process_radtags` with the actual enzyme, barcode file, and barcode
   layout. Inspect its retention report before proceeding.
4. Preserve sample IDs and record input/output read counts. Perform FastQC and
   inspect reports before deciding on adapter and quality filtering.
5. Supply the known library adapter to `workflow.py trim`. This entry point uses
   fastp defaults for quality filtering and an explicit minimum length of 35;
   these are workflow defaults, not a claim about the historical research run.
   Run FastQC again on the resulting reads.

Each trim run writes `reads.fastq.gz` in a sample-specific directory. Before
assembly, stage them under unique sample filenames, for example:

```bash
mkdir -p data/sorted
cp results/trimmed/SAMPLE_A/reads.fastq.gz data/sorted/SAMPLE_A.fastq.gz
```

Never pool files with identical basenames or apply generic trimming to raw
barcode/oligo-bearing reads without checking which bases the protocol needs.

References: [Stacks clone_filter](https://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php),
[Stacks process_radtags](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php),
[fastp documentation](https://github.com/OpenGene/fastp).
