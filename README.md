# Population Structure and Geographic Ancestry of *Ambystoma*

This repository contains the workflow used to estimate and visualize genomic ancestry in 125 *Ambystoma barbouri* and *Ambystoma texanum* individuals sampled across Illinois, Indiana, Kentucky, Missouri, Ohio, Tennessee, and West Virginia.

The analysis converts a filtered RADseq VCF into PLINK format, estimates ancestry with ADMIXTURE for K = 1–10, compares cross-validation errors, and produces STRUCTURE-style barplots and geographic ancestry maps. K = 5 had the lowest cross-validation error, while K = 2, K = 3, and K = 5 were retained for biological interpretation and comparison.

## Repository overview

```text
ANCESTRY125/
├── Amby125.vcf                 # Filtered VCF containing 125 individuals
├── Amby125.csv                 # Original sample coordinates
├── Amby125.bed                 # PLINK binary genotype file
├── Amby125.bim                 # PLINK variant information
├── Amby125.fam                 # PLINK sample information
├── Amby125.cv.error            # Cross-validation errors for K = 1–10
├── Amby125.K.Q                 # Individual ancestry proportions for each K
├── Amby125.K.P                 # Estimated allele frequencies for each K
├── logK.out                    # ADMIXTURE output logs
├── coordinates.csv             # Coordinates formatted for mapmixture
├── admixture2.csv              # K = 2 plotting data
├── admixture3.csv              # K = 3 plotting data, if generated
├── admixture5.csv              # K = 5 plotting data
└── figures/                    # STRUCTURE barplots and ancestry maps
```

`K` in the filenames above represents the tested number of ancestry clusters.

## Software

- PLINK 1.9
- ADMIXTURE 1.3.0
- R 4.5.3
- R packages: `mapmixture`, `ggplot2`, `gridExtra`, `sf`, `ggspatial`, and `maps`

The analysis was run in a Conda environment:

```bash
conda create -n ancestry_env -c conda-forge -c bioconda \
  plink=1.90b6.21 admixture=1.3.0 r-base r-remotes \
  r-ggplot2 r-gridextra r-sf r-ggspatial r-maps

conda activate ancestry_env
```

`mapmixture` was installed from GitHub after installing its spatial dependencies through Conda:

```r
remotes::install_github(
  "Tom-Jenkins/mapmixture",
  force = TRUE,
  upgrade = "never",
  dependencies = FALSE
)
```

## 1. Convert the VCF to PLINK format

```bash
plink --vcf Amby125.vcf \
  --make-bed \
  --out Amby125 \
  --const-fid \
  --allow-extra-chr
```

Because the data were assembled de novo and do not use standard chromosome names, the first column of the BIM file was replaced with `0` for ADMIXTURE compatibility:

```bash
awk '{$1="0"; print $0}' Amby125.bim > Amby125.bim.tmp
mv Amby125.bim.tmp Amby125.bim
```

## 2. Run ADMIXTURE for K = 1–10

```bash
FILE="Amby125"

for K in {1..10}; do
  admixture --cv ${FILE}.bed ${K} > log${K}.out
done
```

Extract the cross-validation results:

```bash
grep "CV" log*.out
```

The observed errors were:

| K | CV error |
|---:|---------:|
| 1 | 0.41163 |
| 2 | 0.33425 |
| 3 | 0.32844 |
| 4 | 0.31690 |
| **5** | **0.30980** |
| 6 | 0.31577 |
| 7 | 0.31450 |
| 8 | 0.32664 |
| 9 | 0.34848 |
| 10 | 0.34705 |

K = 5 produced the lowest cross-validation error and was therefore the best-supported model among the tested values. Lower values such as K = 2 and K = 3 remain useful for visualizing broader hierarchical structure.

Save the errors in a two-column file:

```bash
grep "CV" log*.out | \
  awk '{print $3,$4}' | \
  sed -e 's/(//;s/)//;s/://;s/K=//' > Amby125.cv.error

sort -n -k1,1 Amby125.cv.error -o Amby125.cv.error
```

## 3. Prepare ancestry files for R

ADMIXTURE writes ancestry proportions to `Amby125.K.Q`. The row order in each Q file corresponds exactly to the sample order in `Amby125.fam`.

### K = 2

```bash
paste -d',' \
  <(awk 'BEGIN{OFS=","}{print $2,$2}' Amby125.fam) \
  <(awk 'BEGIN{OFS=","}{$1=$1; print}' Amby125.2.Q) | \
awk 'BEGIN{print "Site,Ind,Cluster1,Cluster2"}{print}' \
> admixture2.csv
```

### K = 5

```bash
paste -d',' \
  <(awk 'BEGIN{OFS=","}{print $2,$2}' Amby125.fam) \
  <(awk 'BEGIN{OFS=","}{$1=$1; print}' Amby125.5.Q) | \
awk 'BEGIN{print "Site,Ind,Cluster1,Cluster2,Cluster3,Cluster4,Cluster5"}{print}' \
> admixture5.csv
```

The same pattern can be used for K = 3 or another value by changing the Q filename and cluster headers.

## 4. Prepare the coordinate file

The original coordinate columns were renamed to the format required by `mapmixture`:

```bash
awk -F',' 'BEGIN{OFS=","}
NR==1 {print "Site","Lat","Lon"; next}
{print $1,$2,$3}' Amby125.csv > coordinates.csv
```

Confirm that all sample IDs occur in both files:

```bash
comm -3 \
  <(tail -n +2 admixture5.csv | cut -d',' -f1 | sort) \
  <(tail -n +2 coordinates.csv | cut -d',' -f1 | sort)
```

No output indicates a complete match.

## 5. Create the K = 2 STRUCTURE barplot

Each vertical bar represents one individual, and each color represents the proportion of ancestry assigned to one of the two inferred clusters. The full sample codes are displayed below the bars.

```r
library(mapmixture)
library(ggplot2)

admixture2 <- read.csv("admixture2.csv")

structure_K2 <- structure_plot(
  admixture_df = admixture2,
  type = "structure",
  cluster_cols = c("#2166AC", "#D73027"),
  site_dividers = TRUE,
  divider_width = 0.3,
  labels = "individual",
  flip_axis = FALSE
) +
  theme(
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5,
      face = "bold"
    ),
    legend.position = "top"
  ) +
  labs(
    title = "Individual Ancestry Assignment (K = 2)",
    x = NULL,
    y = "Ancestry proportion"
  )

ggsave(
  "K2_structure_individual_codes_bold.png",
  plot = structure_K2,
  width = 28,
  height = 9,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)
```

## 6. Create the K = 5 STRUCTURE barplot

```r
admixture5 <- read.csv("admixture5.csv")

structure_K5 <- structure_plot(
  admixture_df = admixture5,
  type = "structure",
  cluster_cols = c(
    "forestgreen",
    "#6A3D9A",
    "#ABD9E9",
    "#FF7F00",
    "#FFD92F"
  ),
  site_dividers = TRUE,
  divider_width = 0.3,
  labels = "individual",
  flip_axis = FALSE
) +
  theme(
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5,
      face = "bold"
    ),
    legend.position = "top"
  ) +
  labs(
    title = "Individual Ancestry Assignment (K = 5)",
    x = NULL,
    y = "Ancestry proportion"
  )

ggsave(
  "K5_structure_individual_codes_bold.png",
  plot = structure_K5,
  width = 28,
  height = 9,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)
```

## 7. Geographic ancestry mapping

Samples with identical coordinates should be summarized by locality before mapping so that overlapping individuals do not create multiple pies at the same position.

```r
library(dplyr)

coordinates <- read.csv("coordinates.csv")

admixture5$Locality <- sub("_[0-9]+$", "", admixture5$Site)
coordinates$Locality <- sub("_[0-9]+$", "", coordinates$Site)

ancestry_locality <- admixture5 |>
  group_by(Locality) |>
  summarise(across(starts_with("Cluster"), mean), .groups = "drop")

coords_locality <- coordinates |>
  group_by(Locality) |>
  summarise(Lat = mean(Lat), Lon = mean(Lon), .groups = "drop")

admixture_map <- ancestry_locality |>
  transmute(
    Site = Locality,
    Ind = Locality,
    Cluster1,
    Cluster2,
    Cluster3,
    Cluster4,
    Cluster5
  )

coordinates_map <- coords_locality |>
  transmute(Site = Locality, Lat, Lon)
```

The sampling region spans approximately 90.54–81.47°W and 35.65–40.78°N. A slightly expanded plotting boundary can therefore be used:

```r
map_K5 <- mapmixture(
  admixture_df = admixture_map,
  coords_df = coordinates_map,
  cluster_cols = c(
    "forestgreen",
    "#6A3D9A",
    "#ABD9E9",
    "#FF7F00",
    "#FFD92F"
  ),
  cluster_names = paste("Ancestry", 1:5),
  crs = 4326,
  boundary = c(
    xmin = -92.0,
    xmax = -80.0,
    ymin = 34.5,
    ymax = 42.5
  ),
  pie_size = 0.22
) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "Geographic Distribution of K = 5 Ancestry",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(
  "K5_geographic_ancestry_map.png",
  plot = map_K5,
  width = 14,
  height = 9,
  dpi = 400,
  bg = "white"
)
```

For the final figure, a state-boundary basemap should show Illinois, Indiana, Kentucky, Missouri, Ohio, Tennessee, and West Virginia. The seven-state layer can be constructed with `maps` and converted to an `sf` object:

```r
library(maps)
library(sf)

state_names <- c(
  "illinois", "indiana", "kentucky", "missouri",
  "ohio", "tennessee", "west virginia"
)

states7 <- maps::map(
  "state",
  regions = state_names,
  fill = TRUE,
  plot = FALSE
) |>
  sf::st_as_sf()
```

State boundaries, watershed polygons, locality labels, north arrows, and scale bars may be added during final figure refinement. The same ancestry colors must be used consistently across barplots and maps for each K.

## Interpretation

- K represents a hypothesized number of ancestry clusters, not automatically a number of species.
- A single-color bar indicates predominant assignment to one inferred ancestry cluster.
- A multicolored bar indicates shared ancestry or possible admixture.
- K = 5 minimizes cross-validation error for this dataset.
- K = 2 and K = 3 reveal broader hierarchical divisions that may remain biologically informative.
- Biological interpretation should integrate ancestry results with geography, watersheds, phylogenetic relationships, sampling history, and other population-genomic analyses.

## Quality-control checks

Before plotting, verify that the sample order and row counts agree:

```bash
wc -l Amby125.fam Amby125.2.Q Amby125.3.Q Amby125.5.Q
```

Each file should contain 125 rows. In R:

```r
dim(admixture2)
dim(admixture5)
dim(coordinates)
sum(is.na(coordinates))
```

Expected dimensions are 125 × 4 for K = 2, 125 × 7 for K = 5, and 125 × 3 for the coordinates. The coordinate file should contain no missing values.

## Notes on reproducibility

- Do not reorder a Q file independently of its corresponding FAM file.
- The `Site` and `Ind` columns initially contain the same full sample code so individual labels remain traceable.
- For locality-level maps, remove only the terminal individual number when deriving locality codes; retain the original identifiers in the individual-level data.
- The original VCF and coordinate file should remain unchanged. Derived plotting files are written separately.
- ADMIXTURE cluster numbers and colors are arbitrary. Across replicate runs or different K values, clusters should be aligned before making direct color-based comparisons.
- For a robust final analysis, multiple ADMIXTURE replicates per K should be run with different random seeds, followed by assessment of run convergence and cluster alignment.

## Citation

If this workflow is used in a publication, cite the original software and package publications or repositories for PLINK, ADMIXTURE, R, and `mapmixture`.

## Author

**Pauline Owusu-Ansah**  
Ph.D. Candidate in Biology  
Miami University

