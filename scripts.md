# Figure 1
```R
# Load libraries
library(edgeR)
library(tidyverse)
library(pheatmap)

# Load count dataframe
all_ <- readRDS('rnaseq_rawcount.rds')

# Load metadata which is id, age, and sex
id <- read.csv('rnaseq_sample_table.csv', header=FALSE)
names(id) <- c('id', 'age', 'sex')

# make edgeR function
myedgeR <- function(dataframe, group, design){
   count <- mutate_all(dataframe, function(x) as.numeric(as.character(gsub(",", "", x))))
   count <- as.matrix(count)
   group <- factor(group)
   d <- DGEList(counts = count, group = group)     
   keep <- filterByExpr(d, design=design) # Using a design matrix is more accurate
   d <- d[keep, , keep.lib.sizes=FALSE]
   d <- calcNormFactors(d)
   # estimateDisp is necessary for differential expression analysis
   d <- estimateDisp(d, design)
   return(d)
}

                       
# make group

str_split(names(all_), pattern='_', simplify=TRUE) %>% .[,1] -> group

# ★★★ Change group names here ★★★
group[group == "furyo"] <- "NWV"
group[group == "ryokou"] <- "WV"

age <- id$age
sex <- id$sex

# make design matrix                       
design <- model.matrix(~ group + age + sex)


# normalize and dispersion estimation
norm <- myedgeR(all_, group, design)


# calculate logcpm
all_normalized <- cpm(norm, log=TRUE, prior.count=3)


names_group <- c('COV102',
'COV112',
'COV27',
'COV53',
'COV65',
'COV70',
'COV75',
'COV81',
'COV97',
'COV1',
'COV108',
'COV33',
'COV38',
'COV41',
'COV43',
'COV45',
'COV48',
'COV57',
'COV6',
'COV68',
'COV72',
'COV73',
'COV84',
'COV86',
'COV87',
'COV90',
'COV96',
'H1',
'H2',
'H7',
'H8',
'H9',
'H11',
'H14',
'H15',
'H19',
'H20',
'H21',
'H22',
'H23',
'H24',
'H25',
'H16'
)

                       
# Keep group information before changing column names for annotation
annotation_group <- group
colnames(all_normalized) <- names_group                       

# make dendrogram
rho <- cor(all_normalized, method="spearman")
d <- as.dist(1-rho)
h <- as.dendrogram(hclust(d, method = "ward.D"))

fit <- glmQLFit(norm, design)

# Specify coefficients related to group (2nd and 3rd) from design matrix column names
qlf <- glmQLFTest(fit, coef=2:3)

top_genes <- topTags(qlf, n=200)
all_top_gene_names <- rownames(top_genes$table)
filtered_gene_names <- all_top_gene_names[!grepl("^IG|^TR", all_top_gene_names)]
top_gene_names <- filtered_gene_names[1:50]



# 3. Prepare data for heatmap
heatmap_data <- all_normalized[top_gene_names, ]


# 4. Create sample annotation (group information)
annotation_col <- data.frame(
  Group = factor(annotation_group)
)

rownames(annotation_col) <- colnames(heatmap_data)


# 5. Create heatmap

pheatmap(
  heatmap_data,
  annotation_col = annotation_col, # Column (sample) annotation
  cluster_cols = as.hclust(h),     # Use the original dendrogram as is
  scale = "row",                    # Scale along rows (Z-score) for each gene
  show_rownames = TRUE,             # Display gene names
  show_colnames = TRUE,             # Display sample names
  main = "Top 50 Differentially Expressed Genes"
)
```

# Figure2
```R
# Load libraries (install pheatmap if not already installed)
library(edgeR)
library(dplyr)
library(pheatmap)

id <- readRDS('singlecell_sample_id.rds')
names(id) <- c('id', 'sex', 'age')
data <- readRDS('singlecell_rawcount.rds')

# Function to perform normalization using edgeR
myedgeR <- function (dataframe, group, design)
{
   count <- mutate_all(dataframe, function(x) as.numeric(as.character(gsub(",","", x))))
   count <- as.matrix(count)
   group <- factor(group)
   # Create DGEList object
   d <- DGEList(counts = count, group = group) 
   # Filtering low-expressed genes
   keep <- filterByExpr(d, design = design) # Using a design matrix is more accurate
   d <- d[keep, , keep.lib.sizes = FALSE]
   d <- calcNormFactors(d)
   d <- estimateDisp(d, design)
   return(d)
}


group <- c('healthy', 'healthy', 'mild', 'mild', 'severe', 'severe', 'healthy','healthy','healthy','healthy', 'severe')
age <- id$age
sex <- id$sex

# Create design matrix
design <- model.matrix(~ group + age + sex)
# Perform normalization and dispersion estimation
dge_obj <- myedgeR(data, group, design)

# Get normalized counts (CPM)
normalized <- cpm(dge_obj, log=TRUE, prior.count=3)

rho <- cor(normalized, method = "spearman")
d <- as.dist(1 - rho)
h <- as.dendrogram(hclust(d, method = "ward.D"))

fit <- glmQLFit(dge_obj, design)
qlf <- glmQLFTest(fit, coef=2:3)




top_genes <- topTags(qlf, n=200)
all_top_gene_names <- rownames(top_genes$table)
filtered_gene_names <- all_top_gene_names[!grepl("^IG|^TR", all_top_gene_names)]
top_gene_names <- filtered_gene_names[1:50]


# 3. Prepare data for heatmap
heatmap_data <- normalized[top_gene_names, ]




# 4. Create sample annotations (group information)
annotation_col <- data.frame(
 Group = factor(group)
)
rownames(annotation_col) <- colnames(heatmap_data)




# 5. Create heatmap


p2 <-pheatmap(
 heatmap_data,
 annotation_col = annotation_col, # Column (sample) annotations
 cluster_cols = as.hclust(h),     # Use the original dendrogram as is
 scale = "row",                    # Scale (Z-score) by row for each gene
 show_rownames = TRUE,             # Display gene names
 show_colnames = TRUE,             # Display sample names
 main = "Top 50 Differentially Expressed Genes (IG/TR removed)"
)
```

# Figure3
```R
# Download the seurat object below. I cannot upload because this file is big.
https://drive.google.com/file/d/1mydb6F2V_z8HRQMxFTJFi1kmNJ4QqfXO/view?usp=sharing
data <- readRDS('240117_jamb_sctransform.rds')

library(Seurat)
data <- subset(data, subset = !is.na(sample_named))
markers <- FindMarkers(
   x,
   ident.1 = "critical",
   ident.2 = "severe",
   group.by = 'condition1',
   test.use = "MAST", 
   latent.vars = 'sample_named')


library(clusterProfiler)
source('Rscript/rprofile.r')
genes <- marker %>% dplyr::filter(p_val_adj < 0.1 & avg_log2FC > 0) %>% rownames()
source("Rscript/convert_geneID.r")
convert_geneID(genes, "SYMBOL", "ENTREZID") -> entre
entre$ENTREZID -> entre
ego2 <- enrichGO(gene         = entre,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
barplot(ego2)
```
# Figure4
```R
# Download the seurat object below. I cannot upload because this file is big.
https://drive.google.com/file/d/1mydb6F2V_z8HRQMxFTJFi1kmNJ4QqfXO/view?usp=sharing
data <- readRDS('240117_jamb_sctransform.rds')
Idents(data) <- 'condition4'
SplitObject(data, split.by='condition4') -> data2
markerlist <- c()
map(data2, function(x){try(
   FindMarkers(
       x,
       ident.1 = "critical",
       ident.2 = "severe",
       group.by = 'condition1',
       test.use = "MAST",
       latent.vars = 'sample_named'))}) -> markerlist

d <- function(x){
   x <- tibble::rownames_to_column(x, 'gene')
   x <- dplyr::filter(x, !grepl('IG', gene))
   return(x)
}

map(markerlist, d) -> markerlist
lists <- c(4332,2357,6282,567,25801,6279,7305,2495,2512,3579,2215,6280,6283,978,3107,3106,5265,1520,23406,9535,2207,6402,25798,2212,6386,728,8826,5879,391,976,11031,10092,7097,6813,5788,11025,29108,10487,1604,1535,1992,387,527,3689,8635,226,1265)


source('scripts/convert_geneID.r')
convert_geneID(lists, 'ENTREZID', 'SYMBOL')$SYMBOL -> neutrophil_related_genes
df <- data.frame()
df <- as.data.frame(neutrophil_related_genes)
modifym <- function(x){
   mutate(x, foldchange = case_when(
       p_val_adj >= 0.1　~ 0,
       p_val_adj < 0.1 ~ avg_log2FC
   ))    -> x
   x %>% dplyr::filter(gene %in% neutrophil_related_genes) -> x
   x <- x %>% dplyr::select(gene, foldchange)
   return(x)
}




map(markerlist, function(x){
   mutate(x, foldchange = case_when(
       p_val_adj >= 0.1　~ 0,
       p_val_adj < 0.1 ~ avg_log2FC
   ))    -> x
   x %>% dplyr::filter(gene %in% neutrophil_related_genes) -> x
   x %>% dplyr::select(gene, foldchange)
}) -> marker_modify2


dplyr::left_join(df, marker_modify2$'Pre-neu', by=c('neutrophil_related_genes' = 'gene')) -> df
dplyr::left_join(df, marker_modify2$'Pro-neu', by=c('neutrophil_related_genes' = 'gene')) -> df
dplyr::left_join(df, marker_modify2$'Mature', by=c('neutrophil_related_genes' = 'gene')) -> df
names(df) <- c('neutrophil_related_genes', 'Pre-neu', 'Pro-neu', 'Mature')
df <- mutate_all(df, ~replace(., is.na(.), 0))
write.csv(df, 'df_for_heatmap.csv', row.names = FALSE)
```
```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('df_for_heatmap.csv', index_col=1)
df = df.drop('Unnamed: 0', axis=1)
sns.heatmap(data=df, vmin=0, cmap="cividis", annot=True, fmt='.2g')
```



