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

# --- Dendrogram creation (User's original code) ---


rho <- cor(normalized, method = "spearman")
d <- as.dist(1 - rho)
h <- as.dendrogram(hclust(d, method = "ward.D"))

fit <- glmQLFit(dge_obj, design)
qlf <- glmQLFTest(fit, coef=2:3)




top_genes <- topTags(qlf, n=200)
all_top_gene_names <- rownames(top_genes$table)

filtered_gene_names <- all_top_gene_names[!grepl("^IG|^TR", all_top_gene_names)]




# Select top 50 genes from the filtered list


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




