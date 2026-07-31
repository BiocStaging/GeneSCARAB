GeneSCARAB
================

### Authors: Marcos Ramos-González, Emma Serrano-Pérez , Víctor Ramos-González, Mercedes García-González, Francisco J. Romero-Campero

## Introduction

To date, although many recent studies have focused on proposing best
practices for circular analysis (Humphreys & Ruxton, 2017; Landler et
al., 2018, 2020, 2021), their adoption in omics-based chronobiology has
been scarce. Often, the goal of these studies is to determine at what
time of day each biological process occurs, analyzing large volumes of
data. The simplest analysis typically involves estimating the time of
maximum expression (or abundance) in the cycle for each element,
referred as phase of the element. This allows for the grouping of each
element in bins based on discrete intervals for the different phases. In
turn, enrichment analysis using one of the standard bioinformatics
ontologies can be applied to each bin, based on software such as Gene
Set Enrichment Analysis (GSEA) (Subramanian et al., 2005). These assume
monotonic data and do not account for distances in a circular space,
incorrectly measuring those distances and leading to potential errors
when the goal is to determine when a process occurs. Similarly, once the
bins are created, the total number of elements associated with each is
highly variable, so the enrichments may be biased as the background
changes in each bin. This is common, for example, in circadian
transcriptomic studies, in which there exist high and low transcription
time windows, resulting in highly unbalanced bins.

GeneSCARAB (Gene Set Circular Analysis and Rhythms Ascertainment in
Bioprocesses) is designed for analyzing the cohesion and clustering
toward specific points in the cycle of biological sets (genes,
metabolites, …) for chronobiology studies following the proposed best
practices in circular statistics. It employs rigorous tests to ensure
maximum statistical power in the analyses and functions to accurately
plot the results. Also, GeneSCARAB seamlessly integrates into common
pipelines for this kind of studies, using as input standard outputs from
rhythm detection software, and standard Bioconductor annotation packages
to adapt the analysis to any species of interest.

## Installation

GeneSCARAB can be installed as follows:

``` r

#This package is currently under review for inclusion in Bioconductor.

#You can install the development version from GitHub (ramosgonzmarc/GeneSCARAB)
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("ramosgonzmarc/GeneSCARAB")
```

## Quick start

The annotation packages used in this tutorial are available at
*[GeneSCARAB-annot-packages](https://github.com/ramosgonzmarc/GeneSCARAB-annot-packages)*
and can be installed as follows once download and decompressed.

``` r
install.packages("./org.Knitens.eg.db/", repos=NULL, type = "source")
install.packages("./org.Otauri.eg.db/", repos=NULL, type = "source")

```

The first step of this tutorial consist on loading the package and the
associated data.

``` r
library(GeneSCARAB)

data("circa_table_genescarab")
total_phases_table_ld <- data.frame(names = rownames(
    circa_table_genescarab
), phase = as.numeric(
    circa_table_genescarab[["ld.peak.time.hours"]]
))
total_phases_table_sd <- data.frame(names = rownames(
    circa_table_genescarab
), phase = as.numeric(
    circa_table_genescarab[["sd.peak.time.hours"]]
))
```

The basic pipeline implemented in GeneSCARAB consists of taking an
annotation package, a table with the estimated phases for each gene, and
a set of GO terms of interest, and the software will determine the
parameters of the circular distributions for each process, as well as
significance measures to identify rhythmic processes with respect to one
or more diel time points. Circular uniformity testing is performed via
Rayleigh (Rayleigh, 1880), Kuiper (Kuiper, 1960), Hermans-Rasson
(HERMANS & RASSON, 1985; Landler et al., 2019b) and Rao’s spacing (Rao,
1976; Landler et al., 2019a) tests.

For this, we will use data from the paper “Multiomics integration
unveils photoperiodic plasticity in the molecular rhythms of marine
phytoplankton” https://doi.org/10.1093/plcell/koaf033.

To get general information about the package, run:

``` r

?`GeneSCARAB-package`

```

We will first generate a subset of GOs using the annotation package.

``` r
library(org.Otauri.eg.db)
functional_data <- select(org.Otauri.eg.db,
    keys = keys(org.Otauri.eg.db, keytype = "GID"),
    columns = c("GID", "GO")
)
complete_gos <- unique(functional_data$GO)
complete_gos <- complete_gos[!is.na(complete_gos)]

set.seed(2345)
subset_gos <- complete_gos[sample(
    1:length(complete_gos), 50,
    replace = FALSE
)]
```

Then we generate a GO terms list with genes per GO term, including
ancestors of selected GO terms.

``` r
go.list.test <- create_gene_list_go(
    go_vector = subset_gos,
    org.package = "org.Otauri.eg.db",
    go_column = "GO", id_column = "GID"
)
```

We can use this list along with the total phase table to create a list
of gene phases per GO for the short-day (SD) condition, and clean it by
removing GOs with no associated rhythmic genes.

``` r
phases.list.sd <- gene_list_to_phases(
    go.list.test,
    total_phases_table_sd
)
phases.list.sd.clean <- phases.list.sd[which(
    sapply(phases.list.sd, nrow) != 0
)]
```

Based on this list, the `complete_circular_table` function allows the
user to calculate a table showing the circular distribution of the
phases of the elements in each set (genes by GO term in this case) in
terms of quartiles, mean, variance, and magnitude (rho), as well as the
simple and BH-adjusted p-values of the Rayleigh, Kuiper, Hermans-Rasson,
and Rao’s spacing tests against circular uniformity.

This function accepts the list of input phases as key arguments, as well
as three types of additional arguments. The function automatically runs
only the Rayleigh and Kuiper tests to save computation time in
situations where HR and Rao are not necessary. Thus, `force.hr.th`
indicates the p-value threshold reached by Rayleigh or Kuiper below
which the user does not consider it necessary to run HR (because
uniformity would have been already rejected), while
`hr.on.large.sets.th` refers to the maximum size of the sets for which
we want to apply HR. GeneSCARAB runs the tests sequentially, so these
same arguments are used to control the execution of the Rao test, the
last of the three. HR and Rao require permutations, so `iter.hr`and
`iter.rao` controls the number of iterations for each run.

``` r
go.circa.table.sd <- complete_circular_table(
    phase.list = phases.list.sd.clean,
    hr.on.large.sets.th = 200,
    force.hr.th = 0.04,
    rao.on.large.sets.th = 200,
    force.rao.th = 0.04,
    iter.rao = 999, iter.hr = 999
)
```

Some software returns discrete phases rather than continuous ones. In
this case, the `complete_circular_table_grouped` function works the same
as `complete_circular_table`, with the difference that it incorporates
tiebreaking in the tests used as described in Landler et al., 2020) and
also allows the user to control the execution of a permutation-based
Kuiper test.

Now we can explore the results. It is recommended to filter by n \> 2
(or more), as HR may give significant results for gene sets with a small
number of associated rhythmic genes (even 1). Additionally, since the
tests performed when using HR and Rao are typically only for negative
results in Rayleigh and Kuiper, the expected distribution of p-values
will favor high p-values, which can lead to an increase in false
negatives when BH-adjusting them. Therefore, it may be more informative
to filter by p-value in these cases. In any case, this decision is left
to the user, and tests can be forced to be performed always, by setting
p-value thresholds to 0 and size thresholds to a large number.

``` r
go.circa.table.sd <- subset(go.circa.table.sd, n > 2)

hr_boolean_sd <- go.circa.table.sd$hr_p_value < 0.05
hr_boolean_sd[is.na(hr_boolean_sd)] <- FALSE
rao_boolean_sd <- go.circa.table.sd$rao_p_value < 0.05
rao_boolean_sd[is.na(rao_boolean_sd)] <- FALSE

sd_rhythmic <- rownames(go.circa.table.sd)[
    go.circa.table.sd$rayleigh_p_value_adj < 0.05 |
        go.circa.table.sd$kuiper_p_value_adj < 0.05 |
        hr_boolean_sd | rao_boolean_sd
]
sd_rhythmic <- unique(sd_rhythmic[!is.na(sd_rhythmic)])
length(sd_rhythmic)
#> [1] 133
```

We can plot these results with 3 different visualizations. Let’s plot a
sample of 8 biological processes. First a circular boxplot. This
function exports a png containing the plot. MetBrewer’s palettes are
used for colloring.

``` r
ex_plots <- c(
    "GO:0016833", "GO:0015698", "GO:0016790",
    "GO:0055085", "GO:0019205", "GO:0004553", "GO:0008169"
)
plot_table_7 <- go.circa.table.sd[ex_plots, ]
circular_boxplot(plot_table_7,
    color.palette = "Tam"
)
```

<figure>
<img
src="https://github.com/ramosgonzmarc/GeneSCARAB/blob/main/man/figures/sd_boxplot-1.png"
alt="Circular boxplot showing the phase distribution of genes associated to some example GO terms." />
<figcaption aria-hidden="true">Circular boxplot showing the phase
distribution of genes associated to some example GO terms.</figcaption>
</figure>

GeneSCARAB also offers a circular dotplot visualization.

``` r
plot_phase_list_7 <- phases.list.sd.clean[ex_plots]
circular_dotplot(plot_phase_list_7,
    color.palette = "Tam"
)
```

<figure>
<img
src="https://github.com/ramosgonzmarc/GeneSCARAB/blob/main/man/figures/sd_dotplot-1.png"
alt="Circular dotplot showing the phase distribution of genes associated to some example GO terms." />
<figcaption aria-hidden="true">Circular dotplot showing the phase
distribution of genes associated to some example GO terms.</figcaption>
</figure>

And also circular histograms. `nbins` control the number of breaks.

``` r
circular_histogram(plot_phase_list_7, color.palette = "Tam", nbins = 48)
```

<figure>
<img
src="https://github.com/ramosgonzmarc/GeneSCARAB/blob/main/man/figures/sd_hist-1.png"
alt="Circular histogram showing the phase distribution of genes associated to some example GO terms." />
<figcaption aria-hidden="true">Circular histogram showing the phase
distribution of genes associated to some example GO terms.</figcaption>
</figure>

Sets that reject uniformity based on HR or Rao but not by Rayleigh or
Kuiper are probably multimodal. These kind of GOs may be further
analyzed using `multimodal_analysis` function. Taking advantage of
*[OptCirClust](https://CRAN.R-project.org/package=OptCirClust)* (Debnath
& Song, 2021), it takes a phase table from an individual set, which has
already been determined to be non-uniform, and allows the user to
identify its modality, identifying processes with 24, 12, and 8-hour
periods. It also clusters the elements (genes) into each of the
corresponding modes and recalculates their associated circular
statistics tables.

``` r
multi_table <- go.circa.table.sd[which(
    go.circa.table.sd$hr_p_value < 0.05
), ]
multi_go <- rownames(multi_table)[1]
phase_table_multi_go <- phases.list.sd.clean[[multi_go]]
multimodal_list <- multimodal_analysis(phase_table_multi_go)
```

This GO term presents significant bimodality, with circular medians in
ZT4 and ZT15 and high rhos We can also use GeneSCARAB to plot these
clusters. When creating the list, remember that it must be a named list.

``` r
circular_dotplot(
    phases.list = list(
        cluster1 = multimodal_list$cluster_1,
        cluster2 = multimodal_list$cluster_2
    ),
    color.palette = "Tam", tr.height = 0.15
)
```

<figure>
<img
src="https://github.com/ramosgonzmarc/GeneSCARAB/blob/main/man/figures/multi_dotplot-1.png"
alt="Circular dotplot showing two distinct clusters of genes belonging to a bimodal gene set." />
<figcaption aria-hidden="true">Circular dotplot showing two distinct
clusters of genes belonging to a bimodal gene set.</figcaption>
</figure>

## Other ways of creating sets

#### KEGG

Just as in the previous case for GO terms, GeneSCARAB allows you to
automatically generate lists of genes associated with specific KEGG
pathways based on an annotation package and a set of KO terms to include
(or calculate it for all those appearing in the annotation package)
using the `create_gene_list_kegg` function instead of
`create_gene_list_go`. In addition, it allows you to filter pathways
based on whether they belong to bacteria, archaea, plants, animals,
fungi, or protists.

#### Custom gene sets

GeneSCARAB also allows the user to use custom gene lists, for example,
based on associations found in other databases or previous analyses. All
you need to do is define a list structure in which each element
corresponds to a vector of genes (or elements).

``` r
setA <- unique(functional_data$GID)[1:50]
setB <- unique(functional_data$GID)[100:120]
setC <- unique(functional_data$GID)[400:450]

custom.list.test <- list(setA = setA, setB = setB, setC = setC)
```

And, as in the previous cases, the gene_list_to_phases function creates
the phase table for the rhythmic genes contained in each of the custom
gene sets.

``` r
phases.list.custom <- gene_list_to_phases(
    custom.list.test, total_phases_table_ld
)
```

## Comparison of phase distributions of sets under two conditions

Different biological processes or gene sets may exhibit significantly
different circular distributions under different experimental
conditions. GeneSCARAB allows for a direct statistical comparison using
the `test_two_dist` function. This is useful, for example, for
determining advancement or delay in the diel expression of different
processes for different conditions or genotypes.

For example, we will compare data collected under short-day (SD) and
long-day (LD) conditions. For this, we first generate the circular
statistics table for LD condition.

``` r
phases.list.ld <- gene_list_to_phases(
    go.list.test, total_phases_table_ld
)
phases.list.ld.clean <- phases.list.ld[which(
    sapply(phases.list.ld, nrow) != 0
)]

go.circa.table.ld <- complete_circular_table(
    phase.list = phases.list.ld.clean,
    hr.on.large.sets.th = 200,
    force.hr.th = 0.04,
    rao.on.large.sets.th = 200,
    force.rao.th = 0.04,
    iter.rao = 999, iter.hr = 999
)

go.circa.table.ld <- subset(go.circa.table.ld, n > 2)

hr_boolean_ld <- go.circa.table.ld$hr_p_value < 0.05
hr_boolean_ld[is.na(hr_boolean_ld)] <- FALSE
rao_boolean_ld <- go.circa.table.ld$rao_p_value < 0.05
rao_boolean_ld[is.na(rao_boolean_ld)] <- FALSE

ld_rhythmic <- rownames(go.circa.table.ld)[
    go.circa.table.ld$rayleigh_p_value_adj < 0.05 |
        go.circa.table.ld$kuiper_p_value_adj < 0.05 |
        hr_boolean_ld | rao_boolean_ld
]
ld_rhythmic <- unique(ld_rhythmic[!is.na(ld_rhythmic)])
length(ld_rhythmic)
#> [1] 137
```

After this, we compare gene phase distributions under both conditions
for the sets identified as non-uniform under both LD and SD.
`test_two_dist` renders a table with the comparison set-by-set, offering
a p-value and adjusted p-value for each comparison, based on the result
of MANOVA over the sines and cosines of the phases (Landler et al.,
2021, 2022).

``` r
rhythmic_gos <- intersect(ld_rhythmic, sd_rhythmic)

diff_distributed_gos <- test_two_dist(
    phases.list.ld.clean[rhythmic_gos], phases.list.sd.clean[rhythmic_gos]
)
sum(diff_distributed_gos$p_value_adj < 0.05)
#> [1] 82
```

And we can plot one of these.

``` r
go_circa_plot_diff <- rbind(
    go.circa.table.ld["GO:0008026", ],
    go.circa.table.sd["GO:0008026", ]
)
rownames(go_circa_plot_diff) <- c("LD", "SD")

circular_boxplot(go_circa_plot_diff, color.palette = "Austria")
```

<figure>
<img
src="https://github.com/ramosgonzmarc/GeneSCARAB/blob/main/man/figures/compare_boxplot-1.png"
alt="Circular boxplot showing differences in phase distribution of genes associated to GO:0008026 gene set due to photoperiod of entrainment." />
<figcaption aria-hidden="true">Circular boxplot showing differences in
phase distribution of genes associated to <a href="GO:0008026"
class="uri">GO:0008026</a> gene set due to photoperiod of
entrainment.</figcaption>
</figure>

## Identification of sets that do not follow the experimental phase distribution

GeneSCARAB also compares a set’s phases against the overall phase
distribution of the entire study. Biologically, it relates to sets of
elements that display specific and tightly scheduled temporal placement,
suggesting that precise temporal control is likely functionally
important for these sets. This comparison is performed by the
`test_against_gen_dist` function, which also relies on MANOVA test on
the sine and cosine values of the phases of the elements in the set and
of the entire set phases (Landler et al., 2021, 2022).

``` r
not_experimentally_distributed_gos <- test_against_gen_dist(
    phases.list.ld.clean[rhythmic_gos], total_phases_table_ld
)
sum(not_experimentally_distributed_gos$gen_dist_p_value_adj < 0.05)
#> [1] 53
```

## Contribution of each individual element to the complete phase distribution of the set

To determine the relative importance of each element to its set’s
statistics, GeneSCARAB also measures the individual contribution of each
element to the mean and circular dispersion of the complete group using
leave-one-out and ranks them.

``` r
processes_for_gene_contribution <- rhythmic_gos[100:105]
gene_contributions <- gene_contribution_to_set(
    circa_result = go.circa.table.ld[processes_for_gene_contribution, ],
    phase_list = phases.list.ld.clean
)
```

Finally, we can identify the genes that have the greatest impact on the
mean and variance and visualize these genes with respect to the complete
distribution.

``` r
ranked_mean <- sapply(gene_contributions, function(x) rownames(x)[1])
ranked_mean
#>      GO:1902494      GO:1905368      GO:1905369      GO:0004497      GO:0016705      GO:0032451 
#> "ostta09g01450" "ostta05g02530" "ostta05g02530" "ostta15g02680" "ostta18g01610" "ostta11g01340"

ranked_rho <- sapply(gene_contributions, function(x) {
    names(
        which.min(x[, "Rho_rank"])
    )
})
ranked_rho
#>      GO:1902494      GO:1905368      GO:1905369      GO:0004497      GO:0016705      GO:0032451 
#> "ostta01g03170" "ostta13g02870" "ostta13g02870" "ostta04g03310" "ostta03g00650" "ostta11g01340"

gene_of_interest_mean <- subset(
    phases.list.ld.clean$`GO:1902494`,
    names == ranked_mean[1]
)
circular_dotplot(
    phases.list = list(
        ranked_gene = gene_of_interest_mean,
        complete_go = phases.list.ld.clean$`GO:1902494`
    ),
    color.palette = "Austria", tr.height = 0.15
)
```

<figure>
<img
src="https://github.com/ramosgonzmarc/GeneSCARAB/blob/main/man/figures/ranked_contribution-1.png"
alt="Circular dotplot showing the gene whose removal exhibits the highest impact on the mean of the GO:1902494 gene set." />
<figcaption aria-hidden="true">Circular dotplot showing the gene whose
removal exhibits the highest impact on the mean of the <a
href="GO:1902494" class="uri">GO:1902494</a> gene set.</figcaption>
</figure>


## External code

Helper functions `HermansRasson2T`, `HermansRasson2PGroupedRad`, 
`RaoTestValue`, `RaoTestUngroupedRad`, `RaoPGroupedRad` and 
`KuiperPGroupedRad` are adapted from the corresponding code in 
https://link.springer.com/article/10.1007/s00265-020-02881-6 and 
https://link.springer.com/article/10.1186/s40462-019-0160-x and are 
licensed under a Creative Commons Attribution 4.0 International License, 
a copy of which can be found 
[here](https://creativecommons.org/licenses/by/4.0/).

## References

Debnath, T., & Song, M. (2021). Fast Optimal Circular Clustering and
Applications on Round Genomes. IEEE/ACM Transactions on Computational
Biology and Bioinformatics, 18(6), 2061–2071.
<https://doi.org/10.1109/TCBB.2021.3077573>

HERMANS, M., & RASSON, J. P. (1985). A new Sobolev test for uniformity
on the circle. Biometrika, 72(3), 698–702.
<https://doi.org/10.1093/biomet/72.3.698>

Kuiper, N. H. (1960). Tests concerning random points on a circle.
Nederl. Akad. Wetensch. Proc. Ser. A, 63(1), 38–47

Landler, L., Ruxton, G. D., & Malkemper, E. P. (2018). Circular data in
biology: advice for effectively implementing statistical procedures.
Behavioral Ecology and Sociobiology, 72(8), 128.
<https://doi.org/10.1007/s00265-018-2538-y>

Landler, L., Ruxton, G. D., & Malkemper, E. P. (2019a). Circular
statistics meets practical limitations: a simulation-based Rao’s spacing
test for non-continuous data. Movement Ecology, 7(1), 15.
<https://doi.org/10.1186/s40462-019-0160-x>

Landler, L., Ruxton, G. D., & Malkemper, E. P. (2019b). The
Hermans–Rasson test as a powerful alternative to the Rayleigh test for
circular statistics in biology. BMC Ecology, 19(1), 30.
<https://doi.org/10.1186/s12898-019-0246-8>

Landler, L., Ruxton, G. D., & Malkemper, E. P. (2020). Grouped circular
data in biology: advice for effectively implementing statistical
procedures. Behavioral Ecology and Sociobiology, 74(8), 100.
<https://doi.org/10.1007/s00265-020-02881-6>

Landler, L., Ruxton, G. D., & Malkemper, E. P. (2021). Advice on
comparing two independent samples of circular data in biology.
Scientific Reports, 11(1), 20337.
<https://doi.org/10.1038/s41598-021-99299-5>

Landler, L., Ruxton, G. D., & Malkemper, E. P. (2022). The multivariate
analysis of variance as a powerful approach for circular data. Movement
Ecology, 10(1), 21. <https://doi.org/10.1186/s40462-022-00323-8>

Rao, J. S. (1976). Some tests based on arc-lengths for the circle.
Sankhyā: The Indian Journal of Statistics, Series B, 329–338.

Rayleigh, Lord. (1880). XII. On the resultant of a large number of
vibrations of the same pitch and of arbitrary phase. The London,
Edinburgh, and Dublin Philosophical Magazine and Journal of Science,
10(60), 73–78. <https://doi.org/10.1080/14786448008626893>

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Debian GNU/Linux 12 (bookworm)
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.11.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.11.0  LAPACK version 3.11.0
#> 
#> locale:
#>  [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C               LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8     LC_MONETARY=es_ES.UTF-8   
#>  [6] LC_MESSAGES=es_ES.UTF-8    LC_PAPER=es_ES.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Madrid
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] org.Otauri.eg.db_0.1 AnnotationDbi_1.75.0 IRanges_2.47.2       S4Vectors_0.51.3     Biobase_2.73.1       BiocGenerics_0.59.7 
#> [7] generics_0.1.4       GeneSCARAB_0.99.0   
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6          DBI_1.3.0             rlang_1.2.0           magrittr_2.0.5        clue_0.3-68           GetoptLong_1.1.1     
#>   [7] otel_0.2.0            matrixStats_1.5.0     e1071_1.7-17          compiler_4.6.1        RSQLite_3.53.2        png_0.1-9            
#>  [13] vctrs_0.7.3           gsl_2.1-9             pkgconfig_2.0.3       shape_1.4.6.1         crayon_1.5.3          fastmap_1.2.0        
#>  [19] MetBrewer_0.2.0       XVector_0.53.0        energy_1.7-12         Ckmeans.1d.dp_4.3.5   rmarkdown_2.31        bit_4.6.0            
#>  [25] xfun_0.58             Rfast_2.1.5.2         cachem_1.1.0          rmio_0.4.0            jsonlite_2.0.0        blob_1.3.0           
#>  [31] rangen_0.0.1          cluster_2.1.8.2       parallel_4.6.1        R6_2.6.1              RColorBrewer_1.1-3    boot_1.3-32          
#>  [37] Rcpp_1.1.1-1.1        Seqinfo_1.3.0         iterators_1.0.14      knitr_1.51            base64enc_0.1-6       OptCirClust_0.0.4    
#>  [43] tidyselect_1.2.1      rnaturalearth_1.2.0   rstudioapi_0.19.0     yaml_2.3.12           doParallel_1.0.17     codetools_0.2-20     
#>  [49] tibble_3.3.1          KEGGREST_1.53.0       S7_0.2.2              evaluate_1.0.5        sf_1.1-1              units_1.0-1          
#>  [55] proxy_0.4-29          RcppParallel_5.1.11-2 circlize_0.4.18       Biostrings_2.81.3     pillar_1.11.1         circular_0.5-2       
#>  [61] BiocManager_1.30.27   KernSmooth_2.23-26    foreach_1.5.2         bigassertr_0.2.0      ggplot2_4.0.3         scales_1.4.0         
#>  [67] BiocStyle_2.41.0      class_7.3-23          glue_1.8.1            CircMLE_0.3.0         tools_4.6.1           Directional_7.6      
#>  [73] Rnanoflann_0.0.3      mvtnorm_1.4-1         rgl_1.3.36            cowplot_1.2.0         grid_4.6.1            plotrix_3.8-14       
#>  [79] Rfast2_0.1.5.6        rbibutils_2.4.1       colorspace_2.1-2      flock_0.7             cli_3.6.6             zigg_0.0.2           
#>  [85] bigparallelr_0.3.2    ComplexHeatmap_2.29.0 dplyr_1.2.1           gtable_0.3.6          digest_0.6.39         classInt_0.4-11      
#>  [91] rjson_0.2.23          htmlwidgets_1.6.4     farver_2.1.2          memoise_2.0.1         htmltools_0.5.9       lifecycle_1.0.5      
#>  [97] httr_1.4.8            GlobalOptions_0.1.4   GO.db_3.23.1          bigstatsr_1.6.2       bit64_4.8.2
```
