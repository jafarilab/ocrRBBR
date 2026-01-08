## ocrRBBR
Predict OCR-driven Boolean rules from multi-omics datasets

ocrRBBR is an R package to infer Boolean rules linking chromatin accessibility (ATAC-seq peaks) to gene expression (RNA-seq) in both bulk and single-cell multiomic datasets. The package identifies combinations of OCRs (Open Chromatin Regions) that best predict the expression state of a target gene.

Features
- Supports bulk and single-cell multiome data with paired RNA-seq and ATAC-seq per cell or sample.
- Outputs interpretable Boolean rules with associated performance metrics.
- Enables parallel computing for faster processing of large datasets.

<br>

## Data Requirements and Normalization
- RNA-seq (bulk): Data are quantile-normalized or TPM-normalized to adjust for sequencing depth and gene length.
- ATAC-seq (bulk): Signal intensities are quantile-normalized to adjust for sequencing depth differences across samples.
- Single-cell RNA-seq: LogNormalized with a scale factor of 10,000 (using Seurat).
- Single-cell ATAC-seq: Normalized using the ReadsInTSS method to adjust for cell-specific variations in sequencing depth.
- ocrRBBR is data-efficient and works well even with small sample sizes. Unlike neural networks, which require many parameters, ocrRBBR uses ridge regression with fewer parameters, making it suitable for datasets with limited samples. It has been tested on both bulk (85 cell types) and single-cell (9,834 cells) datasets and performs similarly on smaller single-cell datasets. 
- In ocrRBBR, samples do not need to originate from the same tissue or cell type. When samples are from the same cell type or tissue, ocrRBBR partitions them into more homogeneous groups, each associated with a distinct Boolean rule within the inferred rule set.

# Table of Contents
- [Dependency](#Dependency)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Inference of OCR-Driven Boolean Rules in Bulk Multiome Datasets](#Usage1)
  - [Inference of OCR-Driven Boolean Rules in single-cell Multiome Datasets](#Usage2)

# Dependency
Please ensure that the required libraries from the following list are installed and loaded. 
The glmnet package is necessary for fitting ridge regressions. To enable parallel computing in ocrRBBR, the following packages must also be installed: doParallel, parallel, foreach, and doSNOW.

```R
library(doParallel)
library(parallel)
library(foreach)
library(doSNOW)
library(glmnet)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)
library(dplyr)
library(Matrix)
library(rtracklayer)
library(S4Vectors)
library(stats)
library(utils)
```  

<br>

# Installation
```R
install.packages("devtools")
devtools::install_github("CompBioIPM/ocrRBBR")
```
<br>

# Usage
# I. Inference of OCR-Driven Boolean Rules in Bulk Multiome Datasets
A toy dataset example is provided in Data/. Please see following examples for instructions.

`ocrRBBR_single_cell()` predicts Boolean rules for a given gene based on bulk-level multi-omics datasets (RNA-seq and ATAC-seq).
```R
ocrRBBR_bulk(
  rnaseq_data  = rnaseq_data,      # Matrix of RNA-seq gene expression (genes × samples)
  atacseq_data = atacseq_data,     # Matrix of ATAC-seq peak accessibility (peaks × samples)
  gene_name    = "Rag2",           # Gene for which Boolean rules will be inferred
  peak_ids     = peak_ids,         # Candidate regulatory peaks for this gene
  max_feature  = 3,                # Maximum number of OCRs allowed in a Boolean rule
  slope        = 10,               # Slope for sigmoid activation function
  num_cores    = 8                 # Number of parallel workers for computation
)
```

Parameter Descriptions   
```bash
# Required arguments
# rnaseq_data    A numeric matrix of RNA-seq expression values.
#                Rows correspond to genes, columns correspond to cell types or samples.
#                **Note:** *ocrRBBR* was tested using **quantile-normalized RNA-seq data**, but it should also work equally well on **TPM-normalized RNA-seq datasets**, provided the data is appropriately scaled across samples.
#
# atacseq_data   A numeric matrix of ATAC-seq signal intensities.
#                Rows correspond to peaks, columns correspond to cell types or samples.
#                Column names must match those of rnaseq_data.
#                **Note:** Similar to RNA-seq data, *ocrRBBR* is tested using **quantile-normalized ATAC-seq data** but is expected to work with other normalization methods, as long as the data distributions are comparable across samples.
#
# gene_name	     A character string specifying the gene for which to infer Boolean rules.
# peak_ids	     A vector of peak identifiers corresponding to rows in atacseq_data to be used as candidate regulatory regions for gene_name.
#
# Optional arguments
# max_feature    An integer specifying the maximum number of input features allowed in a Boolean rule. Default is 3.
# slope	         The slope parameter for the sigmoid activation function. Default is 10.
# num_cores	     The number of parallel workers to use for computation. Adjust according to your system. Default is NA (automatic selection).
```
<br>

#### `link_peaks_to_tss()` links ATAC-seq peaks to genes based on a user-defined window (±100kb by default) around the TSS (Transcription Start Site).
```bash
linked_peaks <- link_peaks_to_tss(
  gtf_file = gtf_file,          # Path to the GTF file with gene annotations
  peaks_gr = peaks_gr,          # GRanges object containing ATAC-seq peaks
  gene_list = NA,               # Optional: A list of specific genes to link peaks (default is NA, considering all genes)
  tss_window = NA               # Optional: A custom window size around the TSS (default is ±100kb)
)
```
<br>

## Example using a mouse multiome dataset
#### Step 1. Load data
```R
# Load the ATAC-seq data, RNA-seq data, and peak locations
data(multiome_human_mouse)

# List all objects in the current R environment
ls()
[1] "human_atacseq_data" "human_cell_type"    "human_meta_data"    "human_peaks_gr"     "human_rnaseq_data"  "mouse_atacseq_data" "mouse_peaks_gr"     "mouse_rnaseq_data" 

# Inspect the GRanges object peaks_gr, which contains peak consensus scores across mammalian genomes and associated peak p-values.
head(mouse_peaks_gr)
GRanges object with 6 ranges and 3 metadata columns:
      seqnames          ranges strand |    peakID phastCons_scores mlog10_bestPvalue
         <Rle>       <IRanges>  <Rle> | <integer>        <numeric>         <numeric>
  [1]     chr1 3020761-3020811      * |         1             0.00              0.56
  [2]     chr1 3087201-3087251      * |         2             0.00              0.50
  [3]     chr1 3120084-3120134      * |         3             0.07             10.80
  [4]     chr1 3121460-3121510      * |         4             0.15              3.02
  [5]     chr1 3372762-3372812      * |         5             0.03              1.31
  [6]     chr1 3399192-3399242      * |         6             0.06              2.39
  -------
  seqinfo: 44 sequences from an unspecified genome; no seqlengths

# Inspect the atacseq_data matrix, where rows correspond to peaks and columns correspond to cell types.
# The values represent quantile-normalized ATAC-seq signal intensities.
mouse_atacseq_data[1:5, 1:5]
       LTHSC.34-.BM LTHSC.34+.BM STHSC.150-.BM MPP4.135+.BM proB.CLP.BM
278345         4.37         2.57          0.90         4.84        3.59
278346         0.41         0.10          0.90         0.11        3.21
278352         0.41         0.10          1.85         0.83        1.15
278353         2.36         0.71          1.85         0.11        2.52
278354         0.41         1.64          0.90         0.83        0.90

# Inspect the rnaseq_data matrix, where rows correspond to genes and columns correspond to cell types.
# The values represent quantile-normalized RNA-seq signal intensities.
mouse_rnaseq_data[ , 1:5]
     LTHSC.34-.BM LTHSC.34+.BM STHSC.150-.BM MPP4.135+.BM proB.CLP.BM
Rag2     1.020795      1.02126       2.98536     56.39795    483.3949
Spi1   117.386316    199.55519     362.58834    458.98175    327.2625
```

#### Step 2. Generate a list of peaks located within a ±100kb window around each gene.
```R
gene_name <- "Rag2"

linked_peaks <- link_peaks_to_tss(
  gtf_file = system.file("extdata", "gencode.vM25.annotation.sample.gtf", package = "ocrRBBR"),
  peaks_gr = mouse_peaks_gr,
  gene_list = gene_name,
  tss_window = 100000 # ±100kb
)

linked_peaks
# A tibble: 83 × 7
   peak                     gene_id            gene_name gene_type      transcript_id      peak_id min_distance
   <chr>                    <chr>              <chr>     <chr>          <chr>                <int>        <dbl>
 1 chr2:101525902-101525952 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278345       98790.
 2 chr2:101533903-101533953 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278346       90790.
 3 chr2:101537021-101537071 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278347       87672.
 4 chr2:101543845-101543895 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278348       80848.
 5 chr2:101545763-101545813 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278349       78930.
 6 chr2:101549272-101549322 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278350       75420.
 7 chr2:101550344-101550394 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278351       74348.
 8 chr2:101551169-101551219 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278352       73524.
 9 chr2:101551591-101551641 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278353       73102.
10 chr2:101552534-101552584 ENSMUSG00000032864 Rag2      protein_coding ENSMUST00000044031  278354       72158.
# ℹ 73 more rows
# ℹ Use `print(n = ...)` to see more rows
# 
```

#### Step 3. Extract peak GRanges and filter by conservation + pvalue
```R
# -------------------------------
# Get linked peaks for gene
# -------------------------------
linked_peaks_gene <- linked_peaks[linked_peaks$gene_name == gene_name, ]
if(nrow(linked_peaks_gene) == 0){
  stop("No linked peaks found for gene: ", gene_name)
}

peak_ids <- linked_peaks_gene$peak_id

# -------------------------------
# This step removes ATAC-seq peaks with low signal intensities (based on p-values) or peaks not conserved across the mammalian genome.
# This step helps reduce potential false positive predictions by ocrRBBR and can be omitted if desired.
# -------------------------------
peaks_gr_tmp <- mouse_peaks_gr[mouse_peaks_gr$peakID %in% peak_ids]

if(length(peaks_gr_tmp) == 0){
  stop("No peak GRanges found for gene: ", gene_name)
}

phast_median <- median(peaks_gr_tmp$phastCons_scores, na.rm = TRUE)
pval_median  <- median(peaks_gr_tmp$mlog10_bestPvalue, na.rm = TRUE)

peaks_gr_tmp <- peaks_gr_tmp[
  (peaks_gr_tmp$phastCons_scores > phast_median) &
    (peaks_gr_tmp$mlog10_bestPvalue > pval_median)
]

peak_ids <- peaks_gr_tmp$peakID
head(peak_ids)
[1] 278345 278346 278352 278353 278354 278355
```

#### Step 4. Train the model and output the predicted Boolean regulatory rules.
```R
boolean_rules <- ocrRBBR_bulk(mouse_rnaseq_data, mouse_atacseq_data, gene_name, peak_ids, max_feature = 3, slope = 10, num_cores = NA)
Starting processing for gene: Rag2 ...
All input checks passed.
training process started with  8  computing cores
  |====================| 100%


head(boolean_rules)
                                                                                                    Boolean_Rule   R2     BIC                           Rule_Coefficients
1                              [AND(278352,278381,278384),AND(¬278352,278381,278384),AND(278352,¬278381,278384)] 0.79 -306.88 0.46:0.72:0.34:-2.45:-2.17:-0.4:-0.64:-0.05
2 [AND(278362,278381,278398),AND(¬278362,278381,278398),AND(278362,¬278381,278398),AND(¬278362,¬278381,¬278398)] 0.79 -306.71  0.49:0.75:0.55:-1.6:-2.02:-0.48:-0.67:0.05
3 [AND(278381,278390,278398),AND(278381,¬278390,278398),AND(¬278381,¬278390,278398),AND(¬278381,278390,¬278398)] 0.79 -306.69 0.18:-1.59:0.93:-0.35:0.58:0.13:-1.33:-0.85
4                              [AND(278355,278381,278384),AND(¬278355,278381,278384),AND(278355,¬278381,278384)] 0.79 -305.85 0.54:0.75:0.34:-1.99:-1.88:-0.49:-0.6:-0.08
5                                                                                           [AND(278386,278398)] 0.75 -304.72                       0.71:-0.8:-1.04:-0.28
6                                                                                                       [278384] 0.73 -304.67                                  0.54:-0.54


# Top-ranked rule sets (high R², low BIC) highlight minimal yet robust combinations of OCRs that are most predictive of gene regulation.

# Each row represents a Boolean rule set linking chromatin accessibility (OCRs; indexed by peak IDs) to gene expression.
# Boolean_Rule shows the logical combinations of peaks (using AND and NOT, ¬) that best explain the gene’s expression pattern. Multiple AND clauses indicate alternative regulatory configurations that lead to similar transcriptional outcomes.
# R² quantifies how much variance in gene expression is explained by the rule set; higher values indicate stronger explanatory power.
# BIC (Bayesian Information Criterion) balances model fit and complexity; lower values indicate a better trade-off between accuracy and simplicity.
# Rule_Coefficients provide the signed contribution of each Boolean rule in the fitted model. Positive coefficients indicate activating regulatory effects, while negative coefficients suggest repressive or inhibitory influences.
```

<br>
<br>

# II. Inference of OCR-Driven Boolean Rules in single-cell Multiome Datasets
A toy dataset example is provided in Data/. Please see following examples for instructions.

`ocrRBBR_single_cell()` infers OCR-driven Boolean regulatory rules for a given gene using single-cell paired RNA-seq and ATAC-seq multiome data, integrating chromatin accessibility states in gene-flanking regions with gene expression across cells.
```R
ocrRBBR_single_cell(
  rnaseq_data  = rnaseq_data,      # RNA-seq expression matrix (genes × cells)
  atacseq_data = atacseq_data,     # ATAC-seq accessibility matrix (peaks × cells)
  gene_name    = "Rag2",           # Gene for which Boolean rules are inferred
  peak_ids     = peak_ids,         # Candidate regulatory peaks for the gene
  max_feature  = 3,                # Maximum number of OCRs in a Boolean rule
  slope        = 10,               # Slope of the sigmoid activation function
  num_cores    = 8,                # Number of parallel workers
  ESS          = 500,              # Effective sample size
  meta.data    = meta.data         # Per-cell metadata (QC metrics)
)
```

Parameter Descriptions   
```bash
# rnaseq_data    A numeric matrix of RNA-seq expression values (genes × cells).
#                RNA-seq values are assumed to be normalized using Seurat’s LogNormalize method with a scale factor of 10,000:
#                NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
#
# atacseq_data   A numeric matrix of ATAC-seq signal intensities (peaks × cells).
#                ATAC-seq counts are assumed to be normalized per cell using the ReadsInTSS method, where raw peak counts are divided by the number of Tn5 insertions within transcription start site (TSS) regions for each cell.
#                ReadsInTSS values are typically obtained from ArchR and applied as column-wise scaling factors.
#
# gene_name      A character string specifying the target gene for Boolean rule inference.
#
# peak_ids       A vector of peak identifiers corresponding to rows in atacseq_data, defining candidate regulatory regions for the gene.
#
# meta_data      A matrix containing per-cell metadata, including:
#                nCount_RNA
#                nFeature_RNA
#                Percentage of mitochondrial reads
#
# Optional arguments
# max_feature    Maximum number of OCRs allowed in a Boolean rule. Default is 3.
# slope          Slope parameter of the sigmoid activation function used in the model. Default is 10.
# num_cores      Number of parallel workers for computation. Default is NA (automatic selection).
# ESS            Effective sample size of the single-cell data after accounting for noise and cell-to-cell correlation. Default is NA.
```
<br>

#### `link_peaks_to_tss()` links ATAC-seq peaks to genes based on a user-defined window (±100kb by default) around the TSS (Transcription Start Site).
```bash
linked_peaks <- link_peaks_to_tss(
  gtf_file = gtf_file,          # Path to the GTF file with gene annotations
  peaks_gr = peaks_gr,          # GRanges object containing ATAC-seq peaks
  gene_list = NA,               # Optional: A list of specific genes to link peaks (default is NA, considering all genes)
  tss_window = NA               # Optional: A custom window size around the TSS (default is ±100kb)
)
```
<br>

#### `ESS()` estimates the effective sample size (ESS) of single-cell RNA-seq data by accounting for correlation among cells within the same cell type.
```bash
ess <- ESS(
  rnaseq_data = rnaseq_data,   # Normalized RNA-seq expression matrix (genes × cells)
  cell_type   = cell_type      # Cell-level metadata with cell type annotations
)
```

Parameter Descriptions   
```bash
# rnaseq_data must be a numeric matrix with genes as rows and cells as columns. RNA-seq data should be normalized (e.g., Seurat LogNormalize, scale factor = 10,000).
# cell_type is a data frame with row names matching the columns of rnaseq_data and a required column named cell_type indicating cell identities.
```
<br>

## Example using a human single-cell multiome dataset
#### Step 1. Load data
```R
# Load the RData file containing the ATAC-seq data, RNA-seq data, and peak locations
data(multiome_human_mouse)

# Inspect the GRanges object peaks_gr, which contains peak consensus scores across mammalian genomes and associated peak p-values.
head(human_peaks_gr)
GRanges object with 6 ranges and 1 metadata column:
      seqnames        ranges strand |    peakID
         <Rle>     <IRanges>  <Rle> | <integer>
  [1]     chr1   10109-10357      * |         1
  [2]     chr1 180730-181630      * |         2
  [3]     chr1 191491-191736      * |         3
  [4]     chr1 267816-268196      * |         4
  [5]     chr1 586028-586373      * |         5
  [6]     chr1 629721-630172      * |         6
  -------
  seqinfo: 33 sequences from an unspecified genome; no seqlengths

# Inspect the atacseq_data matrix, where rows correspond to peaks and columns correspond to cell types.
# The values represent normalized ATAC-seq counts per cell, using the ReadsInTSS method.
human_atacseq_data[1:5, 1:5]
5 x 5 sparse Matrix of class "dgCMatrix"
      AAACAGCCAAGGAATC-1 AAACAGCCAATCCCTT-1 AAACAGCCAATGCGCT-1 AAACAGCCACACTAAT-1 AAACAGCCACCAACCG-1
83441                  .                  .                  .                  .                  .
83442                  .                  .                  .                  .                  .
83443                  .                  .                  .                  .                  .
83444                  .                  .                  .                  .                  .
83445                  .                  .                  .                  .                  .

# Inspect the rnaseq_data matrix, where rows correspond to genes and columns correspond to cell types.
# The values represent normalized RNA-seq signal intensities, calculated using Seurat's LogNormalize method with a scale factor of 10,000..
human_rnaseq_data[1:5, 1:5]
5 x 5 sparse Matrix of class "dgCMatrix"
            AAACAGCCAAGGAATC-1 AAACAGCCAATCCCTT-1 AAACAGCCAATGCGCT-1 AAACAGCCACACTAAT-1 AAACAGCCACCAACCG-1
MIR1302-2HG                  .                  .                  .                  .                  .
FAM138A                      .                  .                  .                  .                  .
OR4F5                        .                  .                  .                  .                  .
AL627309.1                   .                  .                  .                  .                  .
AL627309.3                   .                  .                  .                  .                  .

# Inspect meta data
head(human_meta_data)
                   nCount_RNA nFeature_RNA percent.mt
AAACAGCCAAGGAATC-1       8380         3308   7.470167
AAACAGCCAATCCCTT-1       3771         1896  10.527711
AAACAGCCAATGCGCT-1       6876         2904   6.457243
AAACAGCCACACTAAT-1       1733          846  18.003462
AAACAGCCACCAACCG-1       5415         2282   6.500462
AAACAGCCAGGATAAC-1       2759         1353   6.922798

# Inspect cell_type
head(human_cell_type)
                              cell_id cell_type
AAACAGCCAAGGAATC-1 AAACAGCCAAGGAATC-1 CD4 Naive
AAACAGCCAATCCCTT-1 AAACAGCCAATCCCTT-1   CD4 TCM
AAACAGCCAATGCGCT-1 AAACAGCCAATGCGCT-1 CD4 Naive
AAACAGCCACACTAAT-1 AAACAGCCACACTAAT-1 CD8 Naive
AAACAGCCACCAACCG-1 AAACAGCCACCAACCG-1 CD8 Naive
AAACAGCCAGGATAAC-1 AAACAGCCAGGATAAC-1 CD4 Naive
```

#### Step 2. Generate a list of peaks located within a ±250kb window around each gene.
```R
# Generate a list of peaks located within a ±100kb window around each gene.
gene_name <- "CD74"

linked_peaks <- link_peaks_to_tss(
  gtf_file = system.file("extdata", "gencode.v48.basic.annotation.sample.gtf", package = "ocrRBBR"),
  peaks_gr = human_peaks_gr,
  gene_list = gene_name, 
  tss_window = 250000 # ±250kb
)

linked_peaks
# A tibble: 43 × 7
   peak                     gene_id         gene_name gene_type      transcript_id   peak_id min_distance
   <chr>                    <chr>           <chr>     <chr>          <chr>             <int>        <dbl>
 1 chr5:150175660-150175697 ENSG00000019582 CD74      protein_coding ENST00000353334   83441      237071 
 2 chr5:150180579-150181363 ENSG00000019582 CD74      protein_coding ENST00000353334   83442      231778.
 3 chr5:150185269-150186122 ENSG00000019582 CD74      protein_coding ENST00000353334   83443      227054 
 4 chr5:150208861-150209820 ENSG00000019582 CD74      protein_coding ENST00000353334   83444      203409 
 5 chr5:150217881-150218594 ENSG00000019582 CD74      protein_coding ENST00000353334   83445      194512 
 6 chr5:150234576-150235054 ENSG00000019582 CD74      protein_coding ENST00000353334   83446      177934.
 7 chr5:150264137-150264569 ENSG00000019582 CD74      protein_coding ENST00000353334   83447      148396.
 8 chr5:150350443-150351368 ENSG00000019582 CD74      protein_coding ENST00000353334   83448       61844 
 9 chr5:150352283-150352763 ENSG00000019582 CD74      protein_coding ENST00000353334   83449       60226.
10 chr5:150356700-150359453 ENSG00000019582 CD74      protein_coding ENST00000353334   83450       54673 
# ℹ 33 more rows
# ℹ Use `print(n = ...)` to see more rows

linked_peaks_gene <- linked_peaks[linked_peaks$gene_name == gene_name, ]
if(nrow(linked_peaks_gene) == 0){
  stop("No linked peaks found for gene: ", gene_name)
}

peak_ids <- linked_peaks_gene$peak_id
head(peak_ids)
[1] 83441 83442 83443 83444 83445 83446
```

#### Step 3. Estimate Effective Sample Size (ESS) from Single-Cell RNA-seq Data (Optional)
```R
ess_value <- ESS(rnaseq_data = human_rnaseq_data, cell_type = human_cell_type)
Processing: B intermediate
Processing: B memory
Processing: B naive
Processing: CD14 Mono
Processing: CD16 Mono
Processing: CD4 Naive
Processing: CD4 Proliferating
Processing: CD4 TCM
Processing: CD4 TEM
Processing: CD8 Naive
Processing: CD8 TCM
Processing: CD8 TEM
Processing: cDC2
Processing: dnT
Processing: gdT
Processing: HSPC
Processing: MAIT
Processing: NK
Processing: NK Proliferating
Processing: NK_CD56bright
Processing: pDC
Processing: Plasmablast
Processing: Treg

```

#### Step 4. Train the model and output the predicted Boolean regulatory rules.
```R
boolean_rules <- ocrRBBR_single_cell(human_rnaseq_data, human_atacseq_data, gene_name, peak_ids, max_feature = NA, slope = 6, num_cores = NA, ESS = ess_value, human_meta_data)
Starting processing for gene: CD74 ...
All input checks passed.
training process started with  8  computing cores
  |====================| 100%

head(boolean_rules)
                                                                                                                Boolean_Rule   R2    BIC                          Rule_Coefficients
1                                                                                                                    [83456] 0.40 220.13                                 0.39:-0.38
2                                                                                       [AND(83456,83458),AND(83456,¬83458)] 0.42 222.68                      0.55:-0.09:0.17:-0.68
3                                                                                       [AND(83456,83460),AND(83456,¬83460)] 0.42 224.48                        0.4:-0.08:0.3:-0.68
4                                                                                       [AND(83456,83475),AND(83456,¬83475)] 0.41 226.15                       0.3:-0.08:0.49:-0.64
5                                                                                       [AND(83456,83482),AND(83456,¬83482)] 0.41 226.65                      0.62:-0.42:0.17:-0.36
6  [AND(83456,83458,83460),AND(¬83456,83458,83460),AND(83456,¬83458,83460),AND(83456,83458,¬83460),AND(83456,¬83458,¬83460)] 0.43 239.37 0.43:0.38:0.26:0.55:-0.53:-0.58:0.01:-0.83

# Top-ranked rule sets (high R², low BIC) highlight minimal yet robust combinations of OCRs that are most predictive of gene regulation.

# Each row represents a Boolean rule set linking chromatin accessibility (OCRs; indexed by peak IDs) to gene expression.
# Boolean_Rule shows the logical combinations of peaks (using AND and NOT, ¬) that best explain the gene’s expression pattern. Multiple AND clauses indicate alternative regulatory configurations that lead to similar transcriptional outcomes.
# R² quantifies how much variance in gene expression is explained by the rule set; higher values indicate stronger explanatory power.
# BIC (Bayesian Information Criterion) balances model fit and complexity; lower values indicate a better trade-off between accuracy and simplicity.
# Rule_Coefficients provide the signed contribution of each Boolean rule in the fitted model. Positive coefficients indicate activating regulatory effects, while negative coefficients suggest repressive or inhibitory influences.

```



