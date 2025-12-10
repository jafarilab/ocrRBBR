## ocrRBBR
Predict OCR-driven Boolean rules from multi-omics datasets

ocrRBBR is an R package to infer Boolean rules linking chromatin accessibility (ATAC-seq peaks) to gene expression (RNA-seq) in both bulk and single-cell multiomic datasets. The package identifies combinations of OCRs (Open Chromatin Regions) that best predict the expression state of a target gene.

Features

1. Predict Boolean regulatory rules using bulk-level multiome datasets.
2. Support for single-cell multiome datasets.
3. Handles RNA-seq and ATAC-seq data with matching cell/sample identities.
4. Outputs interpretable Boolean rules with associated metrics.
5. Parallel computing support for faster processing of large datasets.

# Table of Contents
- [Dependency](#Dependency)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Inference of OCR-Driven Boolean Rules in Bulk Multiome Datasets](#usage1)
  - [Inference of OCR-Driven Boolean Rules in single-cell Multiome Datasets](#usage2)

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
ocrRBBR infers Boolean rules for chromatin accessibility states based on two types of multiomic datasets:
1. Bulk paired ATAC-seq and RNA-seq datasets
2. Single-cell paired ATAC-seq and RNA-seq datasets

A toy dataset example is provided in example/. Please see following examples for instructions.

## Inference of OCR-Driven Boolean Rules in Bulk Multiome Datasets
<details>
  <summary>Click me</summary>

# Run ocrRBBR_bulk to infer OCR-driven Boolean rules for a gene

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

This function predicts Boolean rules for a given gene based on bulk-level multi-omics datasets (RNA-seq and ATAC-seq).

Parameter Descriptions   
```bash
Parameters
Argument	      Description

Required arguments:
rnaseq_data  	  A numeric matrix of RNA-seq expression values. Rows correspond to genes, and columns correspond to cell types or samples.
atacseq_data	  A numeric matrix of ATAC-seq signal intensities. Rows correspond to peaks, and columns correspond to cell types or samples.
gene_name	      A character string specifying the gene for which to infer Boolean rules.
peak_ids	      A vector of peak identifiers corresponding to rows in atacseq_data to be used as candidate regulatory regions for gene_name.

Optional arguments:
max_feature	    An integer specifying the maximum number of input features allowed in a Boolean rule. The default is 3.
slope	          The slope parameter for the sigmoid activation function. Default is 10.
num_cores	      The number of parallel workers to use for computation. Adjust according to your system. Default is NA (automatic selection).
```

Output Columns Explanation
```bash
The result from the ocrRBBR_bulk function contains a table of Boolean rules, and their corresponding metrics. Each row corresponds to a predicted Boolean rule set for the gene. Below is an explanation of each column in the output:

Column Name    	Description
Boolean_Rule	  The Boolean rule for the gene based on chromatin accessibility states. Example: [OR(AND(278352, 278381, 278384), AND(~278352, 278381, 278384))] represents a rule with peaks 278352, 278381, and 278384 interacting in an OR operation.
R2	            The adjusted R-squared value of the Boolean rule, indicating the fit quality. A value closer to 1 indicates a better fit.
BIC	            The Bayesian Information Criterion (BIC) score for the model. Lower BIC values indicate a better model fit, penalizing for complexity (more features).
Weights	        The weights associated with each conjunction in the Boolean rule set. These weights represent how strongly each conjunction contributes to the overall rule. For example, 0.46:0.72:0.34:-2.45:-2.17 represents the weights of different rules within the Boolean rule.
```

# `link_peaks_to_tss()` links ATAC-seq peaks to genes based on a user-defined window (±100kb by default) around the TSS (Transcription Start Site).
```bash
linked_peaks <- link_peaks_to_tss(
  gtf_file = gtf_file,          # Path to the GTF file with gene annotations
  peaks_gr = peaks_gr,          # GRanges object containing ATAC-seq peaks
  gene_list = NA,               # Optional: A list of specific genes to link peaks (default is NA, considering all genes)
  tss_window = NA               # Optional: A custom window size around the TSS (default is ±100kb)
)
```

Parameters
```bash
Argument	Description
gtf_file	Path to the GTF file containing the gene annotations (e.g., from Ensembl or Gencode).
peaks_gr	A GRanges object containing peak regions (ATAC-seq peaks).
gene_list	A character vector of gene names to consider. If left empty (default), all genes will be considered.
tss_window	Peaks located within the default or user-specified window size (in base pairs) around the gene TSS are reported. Default is 100,000 bp.
```

Returns
```bash
Returns ATAC-seq peaks located within the specified window size (in base pairs) around the TSS of each gene.
```

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
  gtf_file = "D:\\PBMC\\mm10\\gencode.vM25.annotation.gtf",
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
  
</details>


