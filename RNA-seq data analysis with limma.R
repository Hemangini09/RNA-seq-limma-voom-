# install required packages
BiocManager::install("limma")
BiocManager::install("edgeR")

# import neccessary libraries
library(limma)
library(edgeR)
library(dplyr)


# load expression data (TCGA expresssion data)
exprs <- get(load("dataFilt.RData"))

# Load clinical data 
clinical_data <- read.csv("Clinical.csv")
clinical_data <- data.frame(clinical_data, row.names = T)


# rownames (sample id) of clinical data matches with colnames (sample id) of expression data
all(rownames(clinical_data) %in% colnames(exprs))
all(rownames(clinical_data) == colnames(exprs))


# create a DGE object (DGEList function from edgeR)
dge <- DGEList(exprs)
dge <- calcNormFactors(dge) 
dim(dge)

# make braf_status as factor 
group = factor(clinical_data$braf_status)
group = relevel(group, ref="WT")

# Filter lowly- expressed genes
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge <- dge[-drop,] 
head(dge$counts)

# Design model matrix
mm <- model.matrix(~0+group)

# Voom transformation
v <- voom(dge, mm, plot = T)

# fitting model
fit <- lmFit(v, mm)
head(coef(fit))

# make a contrast to compare samples between two condition 
contr <- makeContrasts(groupV600E-groupWT, levels = colnames(coef(fit)))

# Fitting constrast
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# calculate log2 fold change and other parameters to get significant genes
top.table <- topTable(tmp,adjust.method ="fdr",n = Inf)

# how many genes are up-regulated 
length(which(top.table$adj.P.Val < 0.05  & top.table$logFC > 1))

# how many genes are down-regulated
length(which(top.table$adj.P.Val < 0.05  & top.table$logFC < -1))

# save top.table file
write.csv(top.table,file="top.table.csv",row.names =T )
save(top.table, file="top.table.rda")
