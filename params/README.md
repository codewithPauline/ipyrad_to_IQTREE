# Assembly parameters

Generate a parameter file using the ipyrad version installed in your environment:

```bash
ipyrad -n project
```

Edit the generated `params-project.txt` before running assembly. At minimum review:

| Field | Required decision |
| --- | --- |
| assembly_name / project_dir | Unique project name and writable output directory |
| sorted_fastq_path | Demultiplexed sample reads; preserve unique sample filenames |
| raw_fastq_path / barcodes_path | Leave empty when using already demultiplexed reads |
| assembly_method | `denovo` for this workflow |
| datatype / restriction_overhang | Match the actual library preparation |
| clust_threshold | Set intentionally; historical notes disagree (0.85 vs 0.90) |
| filter_adapters | Set based on preprocessing; historical notes disagree (0 vs 2) |
| min_samples_locus | Choose for the study's sampling and missingness design |
| output_formats | Include PHYLIP (`p`) and any other required outputs |
| pop_assign_file | Configure only when group-specific filtering is intended |

The historical parameter file is preserved in `docs/legacy/ipyrad_step1.txt`.
It is not a verified final parameter set. Do not infer library chemistry or final
filtering thresholds from it. No new final research parameters have been invented.
