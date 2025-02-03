STAT 540 Analysis Assignment
================
<YOUR NAME HERE>

## Link to my assignment repository

Fill in the link to your assignment repository here: [STAT 540 Analysis
Assignment](insert%20link%20here)

## Instructions

- You only need to fill in the code chunks and text where indicated with
  **Exercise \#**.
- When answering the questions, it is not sufficient to provide only
  code. Add explanatory text to interpret your results throughout.
- Make sure to add comments to your code to explain what you are doing
  in each step.

## Academic Integrity

If you use any generative AI tools to complete the assignment, please
provide a brief explanation of how you used them:

1.  What tool(s) did you use (e.g. chatGPT, Claude, GitHub co-pilot)?

<YOUR WRITTEN ANSWER HERE>

2.  How did the tool(s) fit into your workflow? What was input into the
    tool?

<YOUR WRITTEN ANSWER HERE>

3.  What did you receive as output?

<YOUR WRITTEN ANSWER HERE>

4.  Did you modify or verify the output?

<YOUR WRITTEN ANSWER HERE>

### Background

The dataset used for this assignment has been published by [Shi et
al. in 2023](https://www.cell.com/ajhg/fulltext/S0002-9297(23)00317-8)
in the American Journal of Human Genetics. You’ll recognize this study,
“Heterozygous deletion of the autism-associated gene CHD8 impairs
synaptic function through widespread changes in gene expression and
chromatin compaction” as the paper you were asked to review for the
Paper Critique assignment. The raw RNA-seq read counts, and binned
ATAC-seq counts for this study are available in GEO.

## Packages required

All of the packages you will need are listed below. If you have never
used them before, you will need to install them using the
`BiocManager::install("<PACKAGE_NAME>")` or
`install.packages("<PACKAGE NAME)` command. Add any additional packages
you find you need for the code you use to answer the questions.

``` r
library(biomaRt)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(ggrepel)
# use ggplot theme_bw
theme_set(theme_bw())
# set seed for reproducibility
set.seed(3846)
```

## Part 1 - Accessing RNA-seq data using GEOquery

Here we’ll read in the RNA-seq dataset. The RNA-seq read counts are read
in from a compressed `.txt` file, which was obtained from the [GEO
entry](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236993).
The sample metadata is read in with a single line of code using the
`getGEO()` function from the `GEOquery` package.

``` r
# download RNA-seq file 
if(!file.exists("GSE236993.txt.gz")) {
  gzpath <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236993/suppl/GSE236993_CHD8_RNAseq_ReadCounts.txt.gz"
  download.file(gzpath, destfile = "GSE236993.txt.gz")
}

# Read in the file
rna <- read.table(gzfile("GSE236993.txt.gz"), header = TRUE, sep = "\t")

# get metadata
metadata <- getGEO("GSE236993", getGPL = FALSE)[[1]] %>% pData()
```

    ## Found 1 file(s)

    ## GSE236993_series_matrix.txt.gz

## Part 2 - Check the data and repackage into a `DESeqDataSet` object

Does the data look as expected? Let’s check the RNA-seq counts first.
How many samples are there? How many genes? What are the column names?
Are there any duplicate gene names?

``` r
# Check first few rows and dimensions
head(rna)
```

    ##       Gene RNA01 RNA02 RNA03 RNA04 RNA05 RNA06 RNA07 RNA08 RNA09
    ## 1     A1BG     0     0     3     1     3     0     4     1     1
    ## 2 A1BG-AS1     0     1     1     5     1     1     3     0     1
    ## 3     A1CF     0     0     0     0     1     1     0     0     1
    ## 4      A2M     0     1     1     4     4     0     8     3     1
    ## 5  A2M-AS1     0     0     0     1     0     1     0     0     0
    ## 6    A2ML1     0    15    11    31    18    21    41     4    21

``` r
dim(rna)
```

    ## [1] 23895    10

``` r
# check number of unique genes
length(unique(rna$Gene))
```

    ## [1] 23893

We notice there are some duplicate gene names. Let’s remove these. Note
that an alternative strategy may be to average the measurements across
duplicate genes, supposing they are measurements of, say, two
transcripts of the same gene. However, then we would not necessarily
have integer counts, which would complicate our DESeq analysis.

``` r
# remove duplicate gene names
rna <- rna[!duplicated(rna$Gene),]
dim(rna)
```

    ## [1] 23893    10

Next, let’s check the RNA-seq metadata. There are quite a lot of columns
in it, so we’ll pull out only the column(s) we need for the analysis.
The title column contains the sample ID and the condition (H66 control
or CHD8 isoform A or B).

``` r
# Check first few columns and dimensions
metadata %>% tibble()
```

    ## # A tibble: 9 × 43
    ##   title              geo_accession status submission_date last_update_date type 
    ##   <chr>              <chr>         <chr>  <chr>           <chr>            <chr>
    ## 1 RNA01_H66_neuron_… GSM7592755    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 2 RNA02_H66_neuron_… GSM7592756    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 3 RNA03_H66_neuron_… GSM7592757    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 4 RNA04_CHD8A_neuro… GSM7592758    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 5 RNA05_CHD8A_neuro… GSM7592759    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 6 RNA06_CHD8A_neuro… GSM7592760    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 7 RNA07_CHD8B_neuro… GSM7592761    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 8 RNA08_CHD8B_neuro… GSM7592762    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## 9 RNA09_CHD8B_neuro… GSM7592763    Publi… Jul 10 2023     Jul 12 2023      SRA  
    ## # ℹ 37 more variables: channel_count <chr>, source_name_ch1 <chr>,
    ## #   organism_ch1 <chr>, characteristics_ch1 <chr>, characteristics_ch1.1 <chr>,
    ## #   characteristics_ch1.2 <chr>, characteristics_ch1.3 <chr>,
    ## #   treatment_protocol_ch1 <chr>, growth_protocol_ch1 <chr>,
    ## #   molecule_ch1 <chr>, extract_protocol_ch1 <chr>,
    ## #   extract_protocol_ch1.1 <chr>, taxid_ch1 <chr>, data_processing <chr>,
    ## #   data_processing.1 <chr>, data_processing.2 <chr>, platform_id <chr>, …

``` r
dim(metadata)
```

    ## [1]  9 43

``` r
# Extract only the columns we need: 
# title - contains the column name from the count matrix, and the condition (H66 control or CHD8 isoform A or B).
# Extract these pieces of information into separate columns
# Rename the H66 Condition to WT (consistent with presentation in the paper)
coldata <- metadata %>%
  select(title) %>%
  separate(title, 
           into = c("ID", "Condition", "Cell Type", "Replicate"), 
           sep = "_", remove = FALSE) %>%
  select(-`Cell Type`) %>%
  mutate(Condition = ifelse(Condition == "H66", "WT", Condition))
coldata
```

    ##                               title    ID Condition Replicate
    ## GSM7592755   RNA01_H66_neuron_brep1 RNA01        WT     brep1
    ## GSM7592756   RNA02_H66_neuron_brep2 RNA02        WT     brep2
    ## GSM7592757   RNA03_H66_neuron_brep3 RNA03        WT     brep3
    ## GSM7592758 RNA04_CHD8A_neuron_brep1 RNA04     CHD8A     brep1
    ## GSM7592759 RNA05_CHD8A_neuron_brep2 RNA05     CHD8A     brep2
    ## GSM7592760 RNA06_CHD8A_neuron_brep3 RNA06     CHD8A     brep3
    ## GSM7592761 RNA07_CHD8B_neuron_brep1 RNA07     CHD8B     brep1
    ## GSM7592762 RNA08_CHD8B_neuron_brep2 RNA08     CHD8B     brep2
    ## GSM7592763 RNA09_CHD8B_neuron_brep3 RNA09     CHD8B     brep3

Here we will repackage the RNA-seq data into a `DESeqDataSet` object.
This object is used in the `DESeq2` package for differential expression
analysis. It is also a handy format that will join together the counts
and the sample metadata. To do so, we first need to match each sample
from the metadata to the columns in count.

``` r
# move gene ids to rownames to convert rna to numeric matrix
rownames(rna) <- rna$Gene
rna <- as.matrix(rna[, -1])

# make sure the column names match the sample names and are in the same order
identical(colnames(rna), coldata$ID)
```

    ## [1] TRUE

``` r
# create a more informative ID that indicates the condition and replicate number
coldata$ID <- paste0(coldata$Condition, "_", coldata$Replicate)
rownames(coldata) <- coldata$ID
colnames(rna) <- coldata$ID
```

Now, we will use the `DESeqDataSetFromMatrix` function to create a
`DESeqDataSet` object. This object will join the count data and the
sample metadata, as well as specify the design formula for the
differential expression analysis. Note that here we will use “~ 0 +
Condition” in the design formula (the cell means parameterization). This
is because we would like to be able to test for differences between
various condition combinations in our analysis (which will use
contrasts).

``` r
# create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = coldata,
                              design = ~ 0 + Condition)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

## Part 3 - Exploratory data analysis (EDA)

Now that we have the data in the right format, let’s do some exploratory
data analysis.

Recall from the published paper that the authors stated that they
analyzed 15,574 genes, but here we have 23893 genes. Let’s first try to
remove any genes with no expression across samples, and see if it
matches their total.

``` r
# remove genes with no expression
dds <- dds[rowSums(counts(dds)) > 0,]
nrow(dds)
```

    ## [1] 16678

OK, a little closer. Perhaps they used a filter like at least 2 samples
with a count of 1?

``` r
# remove genes with mean count less than 1
dds <- dds[rowSums(counts(dds)) >= 2,]
nrow(dds)
```

    ## [1] 15434

We don’t know for sure what their criteria for “undetected” was since it
wasn’t reported, but this is reasonably close so we’ll move on with this
set of genes.

We’ll next look at the distribution of counts across samples.

``` r
# plot the distribution of counts
assays(dds)$counts %>%
  as.data.frame() %>%
  gather() %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~key) +
  xlab("RNA-seq Counts") +
  ylab("Frequency") 
```

![](analysis_assignment_files/figure-gfm/EDA-1.png)<!-- -->

We can see these distributions are highly skewed, as expected for raw
count data. Let’s remake the previous plot using variance stabilized
counts, which also rescales values to the middle range of sequencing
depth across samples.

``` r
assays(dds)$counts %>%
  vst() %>%
  as.data.frame() %>%
  gather() %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 50) +
  facet_wrap(~key) +
  theme_minimal() +
  xlab("VST normalized RNA-seq Counts") +
  ylab("Frequency") 
```

![](analysis_assignment_files/figure-gfm/EDA%20VST-1.png)<!-- -->

Now For the first Exercise: let’s visualize a sample-sample correlation
heatmap using the variance stabilized counts.

### Exercise 1: Sample-sample correlation heatmap (5 points)

A. Construct a heatmap where each cell represents the spearman
correlation between two samples. Use the variance stabilized counts for
this analysis, and refer back to Lecture 3 for tips on constructing
heatmaps (use an appropriate colour scale, and order rows/columns
meaningfully). (3 points)

``` r
# <YOUR CODE HERE>
```

B. What do you observe that matches your expectations? Do you see
anything surprising? If so, what could explain what you see? (2 points)

<YOUR ANSWER HERE>

## Part 4 - Differential expression analysis

Now that we have the data in the right format and have done some
exploratory data analysis, we can move on to differential expression
analysis. We’ll use `DESeq2` for this analysis, just as the original
authors did.

We want to see if we can reproduce the main results from the paper. The
authors identified 1,169 differentially expressed genes (DEGs) between
CHD8 mutant samples and WT samples (475 increased and 694 decreased in
CHD8 mutant) at a q-value threshold of 0.10.

### Exercise 2: Differential expression analysis with DESeq2 (5 points)

A. Run a DESeq2 analysis comparing the two mutant groups to WT, and
summarize the results. How many genes are differentially expressed at a
BH adjusted p-value threshold of 0.10? How many are increased and how
many are decreased in CHD8 mutant samples compared to WT samples? How
well do these totals align with the original findings? (3 points)

*Hint: to compare the two mutant groups to WT, you can use a numeric
contrast in the `results()` function.*

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER HERE>

B. Diagnostics: Visualize the distribution of p-values, and the
relationship between the mean of normalized counts and shrunken
dispersion estimates. Comment on what you see. (2 points)

``` r
# <YOUR CODE HERE>
```

### Exercise 3: Reproduce the volcano plot (5 points)

Recreate a plot in the style of Figure 2E of the paper (shown below) to
visualize the result you obtained in Exercise 2. Add labels for the top
20 genes by significance. Comment on any major similiarities or
differences.

![](img/fig2E.jpg)

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER HERE>

### Exercise 4: Compare overlap of significant genes (3 points)

The authors provide a list of DEGs in [Supplementary Table
S1](https://www.cell.com/cms/10.1016/j.ajhg.2023.09.004/attachment/38417a36-e543-42e6-a650-dc967c0ea491/mmc2.txt).
This file is provided in this repository for convenience (`mmc2.txt`).

A. How many genes are in this file? Based on the paper, how many genes
do you expect to be in this file? (0.5 point)

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER HERE>

B. Find the proportion of the DEGs you already identified in Exercise 2
that are listed in this file. (0.5 point)

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER HERE>

C. The supplementary file appears to possibly contain a list of DEGs
with adjusted pvalue \< 0.1 for the comparison of CHD8 **isoform B**
samples to WT (notice the range of values in the column that contains
adjusted p-values for isoform B, and that there are no other adjusted
p-values reported). Carry out this comparison (isoform B vs WT) using
contrasts, and find the proportion of these DEGs that are present in
this file. (1 point)

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER HERE>

D. What could the authors do to improve the reproducibility of their
results? (1 point)

<YOUR WRITTEN ANSWER HERE>

### Exercise 5: Visualize results of top DEGs between CHD8 mutant and WT samples (3 points)

Let’s visualize the top 12 DEGs from the comparison of CHD8 mutant
samples to WT samples (from Exercise 2). Use boxplots and/or individual
points to show the distribution of VST counts for these genes across the
three sample groups.

``` r
# <YOUR CODE HERE>
```

### Exercise 6: Visualize results of top DEGs between CHD8 mutant A and CHD8 mutant B samples (4 points)

A. How many DEGs are there between CHD8 mutant A samples to CHD8 mutant
B samples (hint: you’ll need to first run this contrast)? Does this
align with what you expect? (2 points)

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER>

B. Finally, let’s visualize an MA plot of the DE results comparing the
two mutant isoforms in part A, as well as a histogram of the p-values
for this comparison. Comment on what you see. (2 points)

``` r
# <YOUR CODE HERE>
```

<YOUR WRITTEN ANSWER>

## Part 5 - Next steps

From here, we would go on to further investigate the DEGs between CHD8
mutant and WT, perhaps by looking at the pathways they are involved in,
or the biological functions they are associated with, using gene set
enrichment analyses. We would also perform a multi-omic analysis using
the ATAC-seq data to see if there are any changes in chromatin
accessibility that could explain the changes in gene expression. As
these sort of topics come a bit later in the course, they are not
covered in this assignment.
