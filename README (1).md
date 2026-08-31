# Single-End RADseq Pipeline: ipyrad to Phylogenetic Tree

## Overview

This repository documents a single-end RADseq bioinformatics workflow
used to process *Ambystoma* sequencing data and generate a phylogenetic
dataset for downstream population genomic analyses.

The workflow includes:

Raw FASTQ → PCR clone removal → read quality filtering → ipyrad de novo
assembly → locus/SNP filtering → sample quality assessment → removal of
low-coverage individuals → IQ-TREE maximum-likelihood phylogeny

The final filtered dataset contained **114 individuals**.

## 1. Dataset

The original dataset contained:

-   **117 individuals**
-   Single-end RADseq reads
-   Two focal groups:
    -   *Ambystoma barbouri* (Streamside Salamander)
    -   *Ambystoma texanum* (Small-mouthed Salamander)

Analyses were conducted on Linux.

## 2. Clone Filtering

PCR clones were removed from the single-end reads using the Stacks
`clone_filter` program.

``` bash
clone_filter -f sample.fastq -o clone_filtered/
```

Clone filtering was performed for all samples prior to downstream read
processing.

## 3. Read Quality Filtering

Reads were quality-filtered using `fastp`.

``` bash
conda activate fastp_env
```

The clone-filtered `.fq` files were used as input. The resulting cleaned
FASTQ files were subsequently used for the ipyrad assembly.

## 4. ipyrad Assembly

The RADseq data were assembled using **ipyrad** with a **de novo
assembly strategy**.

``` bash
ipyrad -n AMBYSE
nano params-AMBYSE.txt
ipyrad -p params-AMBYSE.txt -s 1234567
```

## 5. Important ipyrad Parameters

Important settings included:

``` text
assembly_method = denovo
clust_threshold = 0.90
filter_adapters = 2
```

The final parameter file should be retained alongside this README as the
authoritative record of all assembly settings.

## 6. Population Assignment File

Each individual was assigned to either `Streamside` or `Smallmouth`.

The working population-specific minimum-sample setting was:

``` text
# Streamside:2 Smallmouth:2
```

This requires a retained locus to contain data from at least 2
Streamside and 2 Smallmouth individuals. This setting was critical for
successful Step 7 filtering.

## 7. Step 7 Filtering

An initial Step 7 run resulted in:

``` text
No loci passed filters
```

Inspection of the ipyrad statistics showed that loci were generated
upstream but were being removed by the population/minimum-sample
filtering criteria.

After changing the population-specific minimum sample requirement to:

``` text
# Streamside:2 Smallmouth:2
```

Step 7 completed successfully.

The assembly contained approximately:

``` text
618,329 prefiltered loci
293,819 loci retained after filtering
```

## 8. ipyrad Output Files

The final ipyrad output directory included formats such as:

``` text
ipyrad_AMBYSE.vcf
ipyrad_AMBYSE.phy
ipyrad_AMBYSE.nex
ipyrad_AMBYSE.str
ipyrad_AMBYSE.ustr
ipyrad_AMBYSE.geno
ipyrad_AMBYSE.ugeno
ipyrad_AMBYSE.usnps
ipyrad_AMBYSE.snpsmap
ipyrad_AMBYSE.loci
ipyrad_AMBYSE.gphocs
ipyrad_AMBYSE.seqs.hdf5
ipyrad_AMBYSE.snps.hdf5
ipyrad_AMBYSE_stats.txt
```

## 9. VCFtools

VCFtools was installed in a dedicated conda environment.

``` bash
conda create -n vcftools_env -c bioconda vcftools
conda activate vcftools_env
vcftools --version
```

## 10. Individual Missingness and Low-Coverage Samples

Individual-level missingness can be examined with:

``` bash
vcftools --vcf ipyrad_AMBYSE.vcf --missing-indv
```

Three individuals showed extremely low locus recovery:

``` text
106     17 loci
301      6 loci
97     145 loci
```

These individuals were excluded from downstream analyses.

``` text
Original samples: 117
Removed samples:    3
Final samples:    114
```

## 11. Removing Low-Coverage Individuals

Create the removal list:

``` bash
echo -e "106\n301\n97" > remove_samples.txt
```

Filter the VCF:

``` bash
vcftools --vcf ipyrad_AMBYSE.vcf \
  --remove remove_samples.txt \
  --recode --recode-INFO-all \
  --out ipyrad_AMBYSE_noLowCov
```

The filtered VCF contains **114 individuals**.

## 12. Filtering Additional ipyrad Outputs

The same individuals were removed from several additional ipyrad
outputs:

``` bash
for f in ipyrad_AMBYSE.ustr ipyrad_AMBYSE.geno ipyrad_AMBYSE.ugeno \
ipyrad_AMBYSE.usnps ipyrad_AMBYSE.snpsmap ipyrad_AMBYSE.loci \
ipyrad_AMBYSE.gphocs; do
grep -vwF -f remove_samples.txt "$f" > "${f%.txt}_noLowCov.txt"
done
```

Filtered datasets use the `ipyrad_AMBYSE_noLowCov` naming convention.

## 13. Filtering the PHYLIP Alignment

Because PHYLIP stores the number of individuals in its header, the
header must also be updated:

``` bash
awk '
NR==FNR {bad[$0]; next}
FNR==1 {sites=$2; next}
{ if (!($1 in bad)) { keep[++n]=$0 } }
END { print n, sites; for (i=1;i<=n;i++) print keep[i] }
' remove_samples.txt ipyrad_AMBYSE.phy > ipyrad_AMBYSE_noLowCov.phy
```

The resulting PHYLIP alignment contains **114 individuals**.

## 14. Extracting Sample IDs from the VCF

Exact sample names and their order in the filtered VCF were extracted
with:

``` bash
bcftools query -l ipyrad_AMBYSE_noLowCov.vcf > sample_list.txt
```

These IDs should be matched exactly to sample metadata used in
downstream analyses.

## 15. IQ-TREE

Phylogenetic reconstruction was performed using **IQ-TREE 3.0.1**.

``` bash
module load iqtree/3.0.1
```

The filtered PHYLIP alignment was used as input:

``` text
ipyrad_AMBYSE_noLowCov.phy
```

## 16. Maximum-Likelihood Tree

IQ-TREE was run with model testing and branch-support analyses:

``` bash
iqtree3 \
  -s ipyrad_AMBYSE_noLowCov.phy \
  -m TEST \
  --alrt 1000 \
  -B 1000 \
  -T AUTO
```

Where:

-   `-m TEST` performs substitution-model testing.
-   `--alrt 1000` performs 1,000 SH-aLRT replicates.
-   `-B 1000` performs 1,000 ultrafast bootstrap replicates.
-   `-T AUTO` detects available CPU threads.

On a multicore Linux workstation or HPC allocation, the thread count can
be specified directly, for example:

``` bash
-T 24
```

## 17. Important IQ-TREE Outputs

Key output files include:

``` text
*.treefile
*.iqtree
*.log
*.ufboot
```

The `*.treefile` contains the inferred maximum-likelihood phylogeny and
can subsequently be visualized and annotated in R using packages such as
`ape`, `ggtree`, and `treeio`.

## 18. Planned Downstream Analyses

The filtered 114-individual dataset can be used for:

-   Maximum-likelihood phylogenetics
-   PCA
-   sNMF / population structure
-   ADMIXTURE
-   Genetic differentiation (FST)
-   Genetic diversity
-   Gene-flow analyses
-   Population assignment
-   Geographic ancestry visualization

## Workflow Summary

``` text
Raw single-end FASTQ
        |
        v
clone_filter
        |
        v
fastp
        |
        v
ipyrad
(de novo assembly)
        |
        v
Steps 1–7
        |
        v
Step 7 locus filtering
        |
        v
293,819 retained loci
        |
        v
VCF / PHYLIP / STRUCTURE / NEXUS outputs
        |
        v
VCFtools QC
        |
        v
Remove low-coverage samples
106, 301, 97
        |
        v
117 → 114 individuals
        |
        v
Filtered PHYLIP
        |
        v
IQ-TREE 3
        |
        v
Maximum-likelihood phylogeny
        |
        v
R visualization + downstream population genomics
```

## Reproducibility

Software used in this workflow includes:

-   Stacks (`clone_filter`)
-   fastp
-   ipyrad
-   VCFtools
-   BCFtools
-   IQ-TREE 3
-   R

Raw sequencing data and large intermediate genomic files are not
included in this repository. Scripts, parameter descriptions, workflow
documentation, and reproducible analysis code can be maintained here
without uploading large FASTQ/VCF datasets.

## Author

**Pauline Owusu-Ansah**\
Ph.D. Biology\
Miami University
