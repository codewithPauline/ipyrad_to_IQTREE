# Single-End RADseq Assembly and Phylogenetic Analysis

## Overview

This repository documents a reproducible **single-end RADseq workflow**
for *Ambystoma* samples, beginning with sequence processing and **de
novo assembly in ipyrad** and ending with **maximum-likelihood
phylogenetic inference in IQ-TREE**.

The overall workflow was:

``` text
Raw single-end RADseq reads
        ↓
Read preprocessing / quality control
        ↓
ipyrad de novo assembly
        ↓
Locus clustering and filtering
        ↓
Final RADseq datasets
        ↓
PHYLIP alignment
        ↓
IQ-TREE
        ↓
Maximum-likelihood phylogeny
```

All samples in the final assembly passed the filtering criteria and were
retained for the phylogenetic analysis. No individuals were removed
between the final ipyrad assembly and IQ-TREE analysis.

------------------------------------------------------------------------

## Study System

The dataset contains RADseq data from salamanders in the genus
*Ambystoma*, including:

-   *Ambystoma barbouri* --- Streamside Salamander
-   *Ambystoma texanum* --- Small-mouthed Salamander

The workflow was designed to generate genomic datasets suitable for
investigating evolutionary relationships and downstream
population-genomic patterns.

------------------------------------------------------------------------

## Software

The principal software used in the workflow included:

-   **Stacks** --- PCR clone filtering
-   **fastp** --- read quality filtering
-   **ipyrad** --- de novo RADseq assembly
-   **IQ-TREE 3** --- maximum-likelihood phylogenetic inference
-   **R** --- downstream visualization and analysis

Conda environments were used where appropriate to maintain reproducible
software installations.

------------------------------------------------------------------------

# 1. Read Preprocessing

Raw single-end RADseq reads were processed before assembly.

PCR duplicates/clones were filtered using the Stacks `clone_filter`
utility.

Example:

``` bash
clone_filter \
    -f sample.fastq \
    -o clone_filtered/
```

The resulting reads were subsequently subjected to quality filtering.

------------------------------------------------------------------------

# 2. Read Quality Filtering

Sequence quality filtering was performed using **fastp**.

Example environment activation:

``` bash
conda activate fastp_env
```

Quality-filtered reads were used as input for the ipyrad assembly.

This preprocessing step was intended to minimize the contribution of
low-quality reads and technical artifacts to locus assembly.

------------------------------------------------------------------------

# 3. ipyrad De Novo Assembly

RADseq loci were assembled using **ipyrad** with a de novo assembly
strategy.

A new ipyrad parameter file was generated using:

``` bash
ipyrad -n AMBYSE
```

The parameter file was then configured for the dataset:

``` bash
nano params-AMBYSE.txt
```

The assembly was run through ipyrad Steps 1--7:

``` bash
ipyrad -p params-AMBYSE.txt -s 1234567
```

The major stages of the ipyrad workflow include:

``` text
Step 1  → Demultiplexing / input organization
Step 2  → Filtering and editing reads
Step 3  → Within-sample clustering
Step 4  → Joint estimation of heterozygosity and sequencing error
Step 5  → Consensus sequence generation
Step 6  → Across-sample clustering
Step 7  → Final locus filtering and output generation
```

------------------------------------------------------------------------

# 4. Assembly Strategy

The assembly used a **de novo** approach rather than mapping reads to a
reference genome.

An important clustering setting was:

``` text
assembly_method = denovo
clust_threshold = 0.90
```

Adapter filtering was also enabled during assembly:

``` text
filter_adapters = 2
```

The complete `params-AMBYSE.txt` file should be retained with the
project to provide the authoritative record of all ipyrad parameters.

------------------------------------------------------------------------

# 5. Population Assignment

Samples were organized into the focal biological groups represented in
the dataset.

Population-specific minimum sampling criteria were incorporated into the
final locus filtering stage.

The working population requirement included:

``` text
Streamside:2
Smallmouth:2
```

This ensured that retained loci contained sufficient representation from
the focal groups for downstream comparative analyses.

------------------------------------------------------------------------

# 6. Final ipyrad Filtering

Step 7 applied the final locus-level filters and generated the
analysis-ready datasets.

The assembly contained approximately:

``` text
Total prefiltered loci: 618,329
Final retained loci:    293,819
```

The filtering procedure removed loci that did not satisfy the specified
assembly criteria while retaining a large genomic dataset for downstream
analysis.

Importantly, **all individuals in the final dataset passed the
sample-level criteria used for this analysis**.

Therefore:

``` text
Final ipyrad samples
        ↓
No sample exclusions
        ↓
All samples retained for phylogenetic inference
```

No post-ipyrad individual removal step was performed for the tree
presented in this workflow.

------------------------------------------------------------------------

# 7. ipyrad Output Files

ipyrad generated multiple output formats that can be used for different
downstream genomic analyses.

Representative outputs include:

``` text
ipyrad_AMBYSE.vcf
ipyrad_AMBYSE.phy
ipyrad_AMBYSE.nex
ipyrad_AMBYSE.str
ipyrad_AMBYSE.ustr
ipyrad_AMBYSE.geno
ipyrad_AMBYSE.ugeno
ipyrad_AMBYSE.snps
ipyrad_AMBYSE.usnps
ipyrad_AMBYSE.snpsmap
ipyrad_AMBYSE.loci
ipyrad_AMBYSE.gphocs
ipyrad_AMBYSE.seqs.hdf5
ipyrad_AMBYSE.snps.hdf5
ipyrad_AMBYSE.stats.txt
```

These formats support multiple downstream applications, including
population structure, genetic differentiation, SNP-based analyses, and
phylogenetics.

------------------------------------------------------------------------

# 8. PHYLIP Alignment for Phylogenetic Analysis

The PHYLIP alignment generated by ipyrad was used for maximum-likelihood
phylogenetic inference.

``` text
ipyrad_AMBYSE.phy
```

Because all final samples were retained, the ipyrad-generated alignment
could be carried directly into the phylogenetic workflow without an
additional sample-removal step.

Before tree inference, the alignment and sample names should be checked
to confirm that:

-   the expected individuals are present;
-   sample identifiers are unique;
-   the PHYLIP header is consistent with the alignment;
-   sample IDs correspond to the project metadata.

------------------------------------------------------------------------

# 9. IQ-TREE Phylogenetic Analysis

Maximum-likelihood phylogenetic inference was performed using **IQ-TREE
3**.

The installation/version can be checked with:

``` bash
iqtree3 --version
```

The ipyrad PHYLIP alignment was supplied directly to IQ-TREE.

A representative analysis was run as:

``` bash
iqtree3 \
    -s ipyrad_AMBYSE.phy \
    -m TEST \
    --alrt 1000 \
    -B 1000 \
    -T AUTO
```

For a system or compute allocation with a known number of available CPU
cores, the number of threads can instead be specified explicitly.

For example:

``` bash
iqtree3 \
    -s ipyrad_AMBYSE.phy \
    -m TEST \
    --alrt 1000 \
    -B 1000 \
    -T 24
```

------------------------------------------------------------------------

# 10. IQ-TREE Options

### Model testing

``` text
-m TEST
```

Tests nucleotide substitution models so that phylogenetic inference can
be performed using an appropriate model for the alignment.

### SH-aLRT branch support

``` text
--alrt 1000
```

Performs **1,000 SH-aLRT replicates** to evaluate branch support.

### Ultrafast bootstrap

``` text
-B 1000
```

Performs **1,000 ultrafast bootstrap replicates**.

### CPU threads

``` text
-T AUTO
```

Allows IQ-TREE to determine the available number of CPU threads.

Alternatively:

``` text
-T 24
```

uses 24 threads when 24 CPU cores have been allocated to the analysis.

------------------------------------------------------------------------

# 11. IQ-TREE Outputs

IQ-TREE generates several files documenting the analysis.

Important outputs include:

``` text
*.treefile
*.iqtree
*.log
*.ufboot
```

### `.treefile`

Contains the final maximum-likelihood tree.

### `.iqtree`

Contains a detailed report of the phylogenetic analysis, including model
information, likelihood statistics, alignment information, and
branch-support results.

### `.log`

Records the progress and commands associated with the IQ-TREE run.

### `.ufboot`

Contains trees generated during the ultrafast bootstrap analysis.

The final `*.treefile` can be imported into R for rooting, annotation,
visualization, and preparation of publication-quality figures.

------------------------------------------------------------------------

# 12. Tree Visualization

The IQ-TREE output can be read into R using packages such as:

``` r
library(ape)
library(ggtree)
library(treeio)
```

For example:

``` r
library(ape)

tree <- read.tree("AMBYSE_iqtree.treefile")

plot(tree, cex = 0.5)
```

Sample metadata can subsequently be matched to tree-tip labels to
visualize geographic locality, population, species assignment, or other
biological variables.

------------------------------------------------------------------------

# 13. Downstream Applications

The ipyrad datasets generated by this workflow can support additional
analyses, including:

-   Maximum-likelihood phylogenetics
-   Population structure analysis
-   sNMF
-   PCA
-   Genetic differentiation (FST)
-   Genetic diversity
-   Admixture analysis
-   Gene-flow analyses
-   Population assignment
-   Geographic ancestry visualization
-   Landscape and population genomics

------------------------------------------------------------------------

# Complete Workflow

``` text
                 SINGLE-END RADseq
                        │
                        ▼
                  Raw FASTQ files
                        │
                        ▼
              PCR clone filtering
                 (clone_filter)
                        │
                        ▼
               Quality filtering
                    (fastp)
                        │
                        ▼
                 ipyrad assembly
                        │
          ┌─────────────┴─────────────┐
          │                           │
          ▼                           ▼
   Within-sample                 Across-sample
     clustering                    clustering
          │                           │
          └─────────────┬─────────────┘
                        ▼
                 Step 7 filtering
                        │
                        ▼
              293,819 retained loci
                        │
                        ▼
              All samples retained
                        │
                        ▼
         Multiple genomic output files
        VCF / PHYLIP / STRUCTURE / etc.
                        │
                        ▼
              PHYLIP alignment
                        │
                        ▼
                   IQ-TREE 3
                        │
          ┌─────────────┼─────────────┐
          ▼             ▼             ▼
      Model test    SH-aLRT 1000   UFBoot 1000
          │             │             │
          └─────────────┼─────────────┘
                        ▼
             Maximum-likelihood tree
                        │
                        ▼
             Visualization in R
```

------------------------------------------------------------------------

# Reproducibility

To improve reproducibility, the repository can include:

``` text
README.md
params-AMBYSE.txt
scripts/
figures/
```

Large sequencing and genomic data files such as FASTQ, VCF, HDF5, and
large alignments do not need to be stored directly in the GitHub
repository.

The repository instead documents the computational workflow, commands,
parameters, and analysis scripts necessary to reproduce the analysis
when the underlying data are available.

------------------------------------------------------------------------

## Author

**Pauline Owusu-Ansah**\
Ph.D. Biology\
Miami University

**Research areas:** Evolutionary biology, population genomics,
computational biology, speciation, hybridization, and gene flow.
