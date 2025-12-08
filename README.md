## ocrRBBR: an R package for the inference of Regression-Based Boolean Rules in multiomic data.
- This repository contains code and tutorials for executing ocrRBBR.
- Sample datasets to execute ocrRBBR are stored within the Data directory.
- The ocrRBBR package supports parallel execution on multi-CPU platforms, enhancing accessibility for real-world Boolean rule inference applications.

#### Step 1. ocrRBBR installation
The ocrRBBR codes are written in R version 4.1.3 and have been tested in both Windows and Linux environments. 

#### Installation
1. Download the compiled package file `RBBR_0.1.0.tar.gz` from this GitHub page.
2. Install the ocrRBBR package by running the following command in R:
   
```R
install.packages("path/to/RBBR_0.1.0.tar.gz", repos = NULL, type = "source")
```
<br>

#### Dependencies  
Please ensure that you have the following packages installed. The glmnet package is required to fit ridge regressions. In order to run ocrRBBR with parallel computing, the packages doParallel, foreach, and doSNOW need to be installed.

```R
install.packages("glmnet")
install.packages("doParallel")  
install.packages("foreach")
install.packages("doSNOW")
```  

<br>

#### Step 2. Prepare input files

```R
#### Load data
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(RBBR)
library(readxl)

# Load the RData file containing the ATAC-seq data, RNA-seq data, and peak locations
load("mouse_dataset.RData")

# List all objects in the current R environment
ls()
[1] "atacseq_data" "peaks_gr"     "rnaseq_data"

# Inspect the GRanges object peaks_gr, which contains peak consensus scores across mammalian genomes and associated peak p-values.
peaks_gr
GRanges object with 512595 ranges and 3 metadata columns:
           seqnames            ranges strand |    peakID phastCons_scores mlog10_bestPvalue
              <Rle>         <IRanges>  <Rle> | <integer>        <numeric>         <numeric>
       [1]     chr1   3020761-3020811      * |         1             0.00              0.56
       [2]     chr1   3087201-3087251      * |         2             0.00              0.50
       [3]     chr1   3120084-3120134      * |         3             0.07             10.80
       [4]     chr1   3121460-3121510      * |         4             0.15              3.02
       [5]     chr1   3372762-3372812      * |         5             0.03              1.31
       ...      ...               ...    ... .       ...              ...               ...
  [512591]     chrY 90812425-90812475      * |    512591                0              3.99
  [512592]     chrY 90812881-90812931      * |    512592                0              3.21
  [512593]     chrY 90813150-90813200      * |    512593                0              0.69
  [512594]     chrY 90813599-90813649      * |    512594                0              0.60
  [512595]     chrY 90828960-90829010      * |    512595                0              1.41
  -------
  seqinfo: 44 sequences from an unspecified genome; no seqlengths

# Inspect the atacseq_data matrix, where rows correspond to peaks and columns correspond to cell types. The values represent quantile-normalized ATAC-seq signal intensities.
atacseq_data[1:5, 1:5]
  LTHSC.34-.BM LTHSC.34+.BM STHSC.150-.BM MPP4.135+.BM proB.CLP.BM
1         0.41         0.71          0.90         0.11        1.94
2         0.41         1.64          0.90         0.83        0.47
3         2.36         0.10          0.90         0.11        0.47
4         0.41         0.10          0.11         0.11        0.79
5         0.41         0.10          0.11         0.11        0.47

# Inspect the rnaseq_data matrix, where rows correspond to genes and columns correspond to cell types. The values represent quantile-normalized RNA-seq signal intensities.
rnaseq_data[1:5, 1:5]
              LTHSC.34-.BM LTHSC.34+.BM STHSC.150-.BM MPP4.135+.BM proB.CLP.BM
0610005C13Rik     1.096732     1.096732       1.02175     1.021812    1.205236
0610007P14Rik   206.053987   246.105317     192.42464   204.298358  189.759175
0610009B22Rik    78.272059    78.837030      68.84475    76.418169  106.085619
0610009L18Rik     8.577159    16.791386      15.51155    16.947354   10.583704
0610009O20Rik   168.645852   157.926022     155.94164   186.261464  162.584556







peaks_gr
GRanges object with 512595 ranges and 3 metadata columns:
           seqnames            ranges strand |    peakID phastCons_scores mlog10_bestPvalue
              <Rle>         <IRanges>  <Rle> | <integer>        <numeric>         <numeric>
       [1]     chr1   3020761-3020811      * |         1             0.00              0.56
       [2]     chr1   3087201-3087251      * |         2             0.00              0.50
       [3]     chr1   3120084-3120134      * |         3             0.07             10.80
       [4]     chr1   3121460-3121510      * |         4             0.15              3.02
       [5]     chr1   3372762-3372812      * |         5             0.03              1.31
       ...      ...               ...    ... .       ...              ...               ...
  [512591]     chrY 90812425-90812475      * |    512591                0              3.99
  [512592]     chrY 90812881-90812931      * |    512592                0              3.21
  [512593]     chrY 90813150-90813200      * |    512593                0              0.69
  [512594]     chrY 90813599-90813649      * |    512594                0              0.60
  [512595]     chrY 90828960-90829010      * |    512595                0              1.41
  -------
  seqinfo: 44 sequences from an unspecified genome; no seqlengths



atacseq <- as.data.frame(read.csv(file = "ImmGenATAC18_AllOCRsInfo.csv", header= TRUE, check.names = FALSE))
rnaseq <- as.data.frame(read.csv(file = "mmc2.csv", header= TRUE, check.names = FALSE))
mmc1 <- as.data.frame(read_excel("mmc1.xlsx", sheet = 1, col_names = TRUE, col_types = "text"))
cell_type_lineage <- mmc1[ ,c(2,4,5)]

#### Extract shared cell types between ATAC-seq and RNA-seq data
cells_types <- intersect(colnames(atacseq), colnames(rnaseq))
```

#### Step 3. Extract ATAC-seq signal intensities and normalize per peak
```R
atacseq_data <- atacseq[   ,(colnames(atacseq) %in% cells_types)]
peak_names <- rownames(atacseq_data)

atacseq_data <- t(atacseq_data)
colnames(atacseq_data) <- peak_names

atacseq_data <- log(1+atacseq_data,10)
atacseq_data_scaled <- atacseq_data
for(j in 1:ncol(atacseq_data)){
  x <- ( atacseq_data[ ,j] - mean(atacseq_data[ ,j]) )/sd(atacseq_data[ ,j])
  atacseq_data_scaled[ ,j]<- 1/(1+exp(-x))
}
```

#### Step 4. Extract ATAC-seq peaks within ±100 kb of the target gene TSS
```R
gene_id   <- "Rag2"

matched_indices <- grep(paste0("(?<!\\w)", gene_id, "(?!\\w)"), atacseq$genes.within.100Kb, perl = TRUE)
atacseq_gene <- atacseq[matched_indices, ]
atacseq_gene <- atacseq_gene[   ,(colnames(atacseq_gene) %in% cells_types)]
```

#### Step 5. Remove ATAC-seq peaks with low signal intensities (based on p-values) or peaks not conserved across the mammalian genome. This step helps reduce potential false positive predictions by ocrRBBR and can be omitted if desired.
```R
peak_info <- atacseq[rownames(atacseq_gene), 1:8]

a <- median(peak_info$mm10.60way.phastCons_scores)
b <- median(peak_info$`_-log10_bestPvalue`)

peak_info <- peak_info[ ((peak_info$mm10.60way.phastCons_scores>=a)&(peak_info$`_-log10_bestPvalue`>=b)), ]
log_atacseq <- atacseq_data_scaled[ ,row.names(peak_info)]
```

#### Step 6. Extract and rescale RNA-seq data for the target gene across blood cell lineages.
```R
rnaseq_gene <- rnaseq[(rnaseq$X %in% gene_id), ]
rnaseq_gene <- rnaseq_gene[ ,(colnames(rnaseq_gene) %in% cells_types)]
rownames(rnaseq_gene) <- gene_id

log_rnaseq <- log(1+rnaseq_gene,10) - min(log(1+rnaseq_gene,10))
log_rnaseq <- log_rnaseq/quantile(as.numeric(unlist(log_rnaseq)), probs = 0.975 , na.rm = TRUE)

data_scaled <- cbind( log_atacseq, t(log_rnaseq) )
data_scaled <- replace(data_scaled, data_scaled>=1, 0.9999)
data_scaled <- replace(data_scaled, data_scaled<=0, 0.0001)

head(data_scaled)
```

#### Step 7. Train the model and output the predicted Boolean regulatory rules.
```R
rbbr           <- rbbr_train(data_scaled, max_feature = min(3,ncol(data_scaled)-1), mode = "1L", slope = 10, penalty = NA, weight_threshold = NA, num_cores = NA)
training process started with  8  computing cores
  |====================| 100%


head(rbbr$boolean_rules_sorted)

                                                                                                        Boolean_Rule                R2       BIC Input_Size Index             Features
1                              [OR(AND(278352,278381,278384),AND(~278352,278381,278384),AND(278352,~278381,278384))] 0.788633167796687 -306.8806          3   706 278352.278381.278384
2 [OR(AND(278362,278381,278398),AND(~278362,278381,278398),AND(278362,~278381,278398),AND(~278362,~278381,~278398))] 0.788220565514429 -306.7148          3  1435 278362.278381.278398
3 [OR(AND(278381,278390,278398),AND(278381,~278390,278398),AND(~278381,~278390,278398),AND(~278381,278390,~278398))] 0.788150615447657 -306.6867          3  2166 278381.278390.278398
4                              [OR(AND(278355,278381,278384),AND(~278355,278381,278384),AND(278355,~278381,278384))] 0.786059559598788 -305.8519          3  1277 278355.278381.278384
5                                                                                               [AND(278386,278398)] 0.746314992241634 -304.7223          2   275        278386.278398
6                                                                                                           [278384] 0.725155753428715 -304.6730          1    16               278384
  Active_Conjunctions                   Weights Layer1, Sub-Rule1
1                   3 0.46:0.72:0.34:-2.45:-2.17:-0.4:-0.64:-0.05
2                   4  0.49:0.75:0.55:-1.6:-2.02:-0.48:-0.67:0.05
3                   4 0.18:-1.59:0.93:-0.35:0.58:0.13:-1.33:-0.85
4                   3 0.54:0.75:0.34:-1.99:-1.88:-0.49:-0.6:-0.08
5                   1                       0.71:-0.8:-1.04:-0.28
6                   1                                  0.54:-0.54
```






























