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
The glmnet package is necessary for fitting ridge regressions. To enable parallel computing in ocrRBBR, the following packages must also be installed: doParallel, foreach, and doSNOW.

```R
library(doParallel)
library(foreach)
library(doSNOW)
library(glmnet)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(ocrRBBR)
library(readxl)
```  

<br>

# Installation
The ocrRBBR codes are written in R version 4.1.3 and have been tested in both Windows and Linux environments.
1. Download the compiled package file `RBBR_0.1.0.tar.gz` from this GitHub page.
2. Install the ocrRBBR package by running the following command in R:
   
```R
install.packages("path/to/RBBR_0.1.0.tar.gz", repos = NULL, type = "source")
```
<br>

# Usage
# I. Inference of OCR-Driven Boolean Rules in Bulk Multiome Datasets
A toy dataset example is provided in example/. Please see following examples for instructions.

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
# max_feature    An integer specifying the maximum number of input features allowed in a Boolean rule. The default is 3.
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
# Load the RData file containing the ATAC-seq data, RNA-seq data, and peak locations
load("mouse_dataset.RData")

# List all objects in the current R environment
ls()
[1] "atacseq_data" "peaks_gr"     "rnaseq_data"

# Inspect the GRanges object peaks_gr, which contains peak consensus scores across mammalian genomes and associated peak p-values.
head(peaks_gr)
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
atacseq_data[1:5, 1:5]
  LTHSC.34-.BM LTHSC.34+.BM STHSC.150-.BM MPP4.135+.BM proB.CLP.BM
1         0.41         0.71          0.90         0.11        1.94
2         0.41         1.64          0.90         0.83        0.47
3         2.36         0.10          0.90         0.11        0.47
4         0.41         0.10          0.11         0.11        0.79
5         0.41         0.10          0.11         0.11        0.47

# Inspect the rnaseq_data matrix, where rows correspond to genes and columns correspond to cell types.
# The values represent quantile-normalized RNA-seq signal intensities.
rnaseq_data[1:5, 1:5]
              LTHSC.34-.BM LTHSC.34+.BM STHSC.150-.BM MPP4.135+.BM proB.CLP.BM
0610005C13Rik     1.096732     1.096732       1.02175     1.021812    1.205236
0610007P14Rik   206.053987   246.105317     192.42464   204.298358  189.759175
0610009B22Rik    78.272059    78.837030      68.84475    76.418169  106.085619
0610009L18Rik     8.577159    16.791386      15.51155    16.947354   10.583704
0610009O20Rik   168.645852   157.926022     155.94164   186.261464  162.584556
```

#### Step 2. Generate a list of peaks located within a ±100kb window around each gene.
```R
gene_name <- "Rag2"

linked_peaks <- link_peaks_to_tss(
  gtf_file = "\path_to\gencode.vM25.annotation.gtf",
  peaks_gr = peaks_gr,
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
peaks_gr_tmp <- peaks_gr[peaks_gr$peakID %in% peak_ids]

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
```

#### Step 4. Train the model and output the predicted Boolean regulatory rules.
```R
res <- ocrRBBR_bulk(rnaseq_data, atacseq_data, gene_name, peak_ids, max_feature = 3, slope = 10, num_cores = NA)

▶️ Starting processing for gene: Rag2 ...
✔️ All input checks passed.
training process started with  8  computing cores
  |====================| 100%


head(res$boolean_rules_sorted)
                                                                                                        Boolean_Rule                R2       BIC Input_Size Index             Features
1                              [OR(AND(278352,278381,278384),AND(~278352,278381,278384),AND(278352,~278381,278384))] 0.788633167796687 -306.8806          3   587 278352.278381.278384
2 [OR(AND(278362,278381,278398),AND(~278362,278381,278398),AND(278362,~278381,278398),AND(~278362,~278381,~278398))] 0.788220565514429 -306.7148          3  1171 278362.278381.278398
3 [OR(AND(278381,278390,278398),AND(278381,~278390,278398),AND(~278381,~278390,278398),AND(~278381,278390,~278398))] 0.788150615447657 -306.6867          3  1673 278381.278390.278398
4                              [OR(AND(278355,278381,278384),AND(~278355,278381,278384),AND(278355,~278381,278384))] 0.786059559598788 -305.8519          3  1047 278355.278381.278384
5                                                                                               [AND(278386,278398)] 0.746314992241634 -304.7223          2   234        278386.278398
6                                                                                                           [278384] 0.725155753428715 -304.6730          1    15               278384
  Active_Conjunctions                   Weights Layer1, Sub-Rule1
1                   3 0.46:0.72:0.34:-2.45:-2.17:-0.4:-0.64:-0.05
2                   4  0.49:0.75:0.55:-1.6:-2.02:-0.48:-0.67:0.05
3                   4 0.18:-1.59:0.93:-0.35:0.58:0.13:-1.33:-0.85
4                   3 0.54:0.75:0.34:-1.99:-1.88:-0.49:-0.6:-0.08
5                   1                       0.71:-0.8:-1.04:-0.28
6                   1                                  0.54:-0.54
```

<br>
<br>

# II. Inference of OCR-Driven Boolean Rules in single-cell Multiome Datasets
A toy dataset example is provided in example/. Please see following examples for instructions.

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
#                ATAC-seq counts are assumed to be normalized per cell using the ReadsInTSS method, where raw peak counts are divided by the number of Tn5 insertions within transcription start site (TSS) regions.
#                ReadsInTSS values are typically obtained from ArchR and applied as column-wise scaling factors.
#
# gene_name      A character string specifying the target gene for Boolean rule inference.
#
# peak_ids       A vector of peak identifiers corresponding to rows in atacseq_data, defining candidate regulatory regions for the gene.
#
# meta.data      A matrix containing per-cell metadata, including:
#                RNA counts per cell 
#                Number of detected features per cell
#                Percentage of mitochondrial reads
#
# Optional arguments
# max_feature    Maximum number of OCRs allowed in a Boolean rule. Default: 3
# slope          Slope parameter of the sigmoid activation function used in the model. Default: 10
# num_cores      Number of parallel workers for computation. Default: NA (automatic selection)
# ESS            Effective sample size of the single-cell data after accounting for noise and cell-to-cell correlation. Default: NA
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
  cell_data   = cell_data      # Cell-level metadata with cell type annotations
)
```

Parameter Descriptions   
```bash
# rnaseq_data must be a numeric matrix with genes as rows and cells as columns. RNA-seq data should be normalized (e.g., Seurat LogNormalize, scale factor = 10,000).
# cell_data must be a data frame whose row names match the column names of rnaseq_data. cell_data must contain a column named predicted_celltype_l2.
```
<br>

## Example using a human single-cell multiome dataset
#### Step 1. Load data
```R
# Load the RData file containing the ATAC-seq data, RNA-seq data, and peak locations
load("human_dataset.RData")

# List all objects in the current R environment
ls()
[1] "atacseq_data" "cell_data"    "meta.data"    "peaks_gr"     "rnaseq_data" 

# Inspect the GRanges object peaks_gr, which contains peak consensus scores across mammalian genomes and associated peak p-values.
head(peaks_gr)
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
# The values represent quantile-normalized ATAC-seq signal intensities.
atacseq_data[1:5, 1:5]
5 x 5 sparse Matrix of class "dgCMatrix"
  AAACAGCCAAGGAATC-1 AAACAGCCAATCCCTT-1 AAACAGCCAATGCGCT-1 AAACAGCCACACTAAT-1 AAACAGCCACCAACCG-1
1                  .                  .                  .                  .                  .
2                  .                  .                  .                  .                  .
3                  .                  .                  .                  .                  .
4                  .                  .                  .                  .                  .
5                  .                  .                  .                  .                  .

# Inspect the rnaseq_data matrix, where rows correspond to genes and columns correspond to cell types.
# The values represent quantile-normalized RNA-seq signal intensities.
rnaseq_data[1:5, 1:5]
5 x 5 sparse Matrix of class "dgCMatrix"
            AAACAGCCAAGGAATC-1 AAACAGCCAATCCCTT-1 AAACAGCCAATGCGCT-1 AAACAGCCACACTAAT-1 AAACAGCCACCAACCG-1
MIR1302-2HG                  .                  .                  .                  .                  .
FAM138A                      .                  .                  .                  .                  .
OR4F5                        .                  .                  .                  .                  .
AL627309.1                   .                  .                  .                  .                  .
AL627309.3                   .                  .                  .                  .                  .

# Inspect metadata
head(meta.data)
                   nCount_RNA nFeature_RNA percent.mt
AAACAGCCAAGGAATC-1       8380         3308   7.470167
AAACAGCCAATCCCTT-1       3771         1896  10.527711
AAACAGCCAATGCGCT-1       6876         2904   6.457243
AAACAGCCACACTAAT-1       1733          846  18.003462
AAACAGCCACCAACCG-1       5415         2282   6.500462
AAACAGCCAGGATAAC-1       2759         1353   6.922798

# Inspect cell_data
head(cell_data)
                              cell_id predicted_celltype_l1 predicted_celltype_l2
AAACAGCCAAGGAATC-1 AAACAGCCAAGGAATC-1                 CD4 T             CD4 Naive
AAACAGCCAATCCCTT-1 AAACAGCCAATCCCTT-1                 CD4 T               CD4 TCM
AAACAGCCAATGCGCT-1 AAACAGCCAATGCGCT-1                 CD4 T             CD4 Naive
AAACAGCCACACTAAT-1 AAACAGCCACACTAAT-1                 CD8 T             CD8 Naive
AAACAGCCACCAACCG-1 AAACAGCCACCAACCG-1                 CD8 T             CD8 Naive
AAACAGCCAGGATAAC-1 AAACAGCCAGGATAAC-1                 CD4 T             CD4 Naive
```

#### Step 2. Generate a list of peaks located within a ±250kb window around each gene.
```R
# Generate a list of peaks located within a ±100kb window around each gene.
gene_name <- "ZEB2"

linked_peaks <- link_peaks_to_tss(
  gtf_file = "\path_to\gencode.v48.basic.annotation.gtf",
  peaks_gr = peaks_gr,
  gene_list = gene_name, 
  tss_window = 250000 # ±250kb
)

linked_peaks
# A tibble: 81 × 7
   peak                     gene_id         gene_name gene_type      transcript_id   peak_id min_distance
   <chr>                    <chr>           <chr>     <chr>          <chr>             <int>        <dbl>
 1 chr2:144183459-144184399 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57884      247004.
 2 chr2:144190904-144191615 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57885      239674 
 3 chr2:144197988-144198604 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57886      232638.
 4 chr2:144216246-144216612 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57887      214504.
 5 chr2:144227723-144227793 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57888      203176.
 6 chr2:144233311-144233665 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57889      197446.
 7 chr2:144234277-144234821 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57890      196384.
 8 chr2:144237655-144238696 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57891      192758 
 9 chr2:144240752-144241365 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57892      189875 
10 chr2:144267141-144267591 ENSG00000169554 ZEB2      protein_coding ENST00000638087   57893      163568.
# ℹ 71 more rows
# ℹ Use `print(n = ...)` to see more rows
```

#### Step 3. Extract peak GRanges
```R
# -------------------------------
# Get linked peaks for gene
# -------------------------------
linked_peaks_gene <- linked_peaks[linked_peaks$gene_name == gene_name, ]
if(nrow(linked_peaks_gene) == 0){
  stop("No linked peaks found for gene: ", gene_name)
}

peak_ids <- linked_peaks_gene$peak_id
```

#### Step 4. Estimate Effective Sample Size (ESS) from Single-Cell RNA-seq Data
```R
ess_value <- ESS(rnaseq_data = rnaseq_data, cell_data = cell_data)
```

#### Step 5. Train the model and output the predicted Boolean regulatory rules.
```R
res <- ocrRBBR_single_cell(rnaseq_data, atacseq_data, gene_name, peak_ids, max_feature = NA, slope = 6, num_cores = NA, ESS = ess_value, meta.data)

▶️ Starting processing for gene: ZEB2 ...
✔️ All input checks passed.
training process started with  8  computing cores
  |====================| 100%

head(res$boolean_rules_sorted)
                                                                                                                                            Boolean_Rule                R2       BIC Input_Size
1                          [OR(AND(57908,57940,57956),AND(~57908,57940,57956),AND(57908,~57940,57956),AND(57908,57940,~57956),AND(~57908,~57940,57956))] 0.512507251846455 -26136.99          3
2                          [OR(AND(57940,57956,57963),AND(~57940,57956,57963),AND(57940,~57956,57963),AND(57940,57956,~57963),AND(~57940,57956,~57963))] 0.496254863988792 -25814.48          3
3                          [OR(AND(57940,57956,57964),AND(~57940,57956,57964),AND(57940,~57956,57964),AND(57940,57956,~57964),AND(~57940,57956,~57964))] 0.495792481048702 -25805.46          3
4                          [OR(AND(57940,57956,57961),AND(~57940,57956,57961),AND(57940,~57956,57961),AND(57940,57956,~57961),AND(~57940,57956,~57961))] 0.494215674451035 -25774.75          3
5 [OR(AND(57921,57940,57956),AND(~57921,57940,57956),AND(57921,~57940,57956),AND(57921,57940,~57956),AND(~57921,~57940,57956),AND(~57921,57940,~57956))] 0.493175035076764 -25754.54          3
6                          [OR(AND(57940,57950,57956),AND(~57940,57950,57956),AND(57940,~57950,57956),AND(57940,57950,~57956),AND(~57940,~57950,57956))] 0.487814267376566 -25651.07          3
  Index          Features Active_Conjunctions                  Weights Layer1, Sub-Rule1
1 57316 57908.57940.57956                   5 0.23:0.47:0.67:0.67:0.34:-0.13:-0.04:-1.74
2 83267 57940.57956.57963                   5   0.2:0.64:0.65:0.49:-0.14:0.41:-0.1:-1.71
3 83268 57940.57956.57964                   5 0.21:0.65:0.66:0.48:-0.22:0.37:-0.13:-1.62
4 83265 57940.57956.57961                   5    0.2:0.63:0.63:0.5:-0.16:0.4:-0.08:-1.68
5 72695 57921.57940.57956                   6      0.17:0.51:0.74:0.55:0.3:0:-0.19:-1.66
6 83197 57940.57950.57956                   5  0.17:0.6:0.53:0.57:0.44:-0.19:-0.04:-1.67
```



