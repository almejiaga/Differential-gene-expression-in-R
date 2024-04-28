library(argparse)
library(DESeq2)
library(ggplot2)
# Create argument parser
parser <- ArgumentParser(description = "Script for differential expression analysis")

# Add arguments
parser$add_argument("-c", "--count_table", dest = "count_table", 
                    help = "Path to the gene expression matrix file")
parser$add_argument("-f", "--factors_table", dest = "factors_table", 
                    help = "Path to the factors table file")
parser$add_argument("-c1", "--condition1", dest = "condition1", 
                    help = "First condition for analysis")
parser$add_argument("-c2", "--condition2", dest = "condition2", 
                    help = "Second condition for analysis")
parser$add_argument("-l", "--label", dest = "label", 
                    help = "Label to be prepended to output files")

# Parse arguments
args <- parser$parse_args()
# Function to construct output file names with label
construct_output_filename <- function(base_name, label) {
  if (is.null(label) || label == "") {
    return(base_name)
  } else {
    return(paste(label, base_name, sep = "_"))
  }
}

# Example usage of the construct_output_filename function
label <- args$label
#define the names of the output
# Modify file names with label
DEGs_output_file <- construct_output_filename("significant_DEGs.csv", label)
volcanoCOVSCA_output_file <- construct_output_filename("volcanoplot.png", label)
volcanoplotCAVSM1_output_file <- construct_output_filename("enhancedvolcano.png", label)

# Reading the gene expression matrix
countTable <- read.table(args$count_table, header = TRUE, row.names = 1)

# Reading the factors table
colData <- read.table(args$factors_table, header = TRUE)

# Filter to include only two groups for analysis
colData2 <- dplyr::filter(colData, Condition %in% c(args$condition1, args$condition2))
colData2$sample <- as.character(colData2$sample)
#filtering count table to exclude genes with less than 32 counts for any treatment
countTable1 = countTable [,6:26]
countTable1 <- countTable1[,sort(colnames(countTable1))]
vector4 <- seq(1:84)
colnames(countTable1) <- vector4[64:84]
countTable1<- countTable1[, colData2$sample]
countTable1[, "max"]= apply(countTable1[, 1:ncol(countTable1)], 1, max)
countTable1=countTable1[countTable1[,ncol(countTable1)]>32,]
countTable1= countTable1 [,-ncol(countTable1)]
#Building deseq2 object with the expression data and the metadata
dds2<-DESeqDataSetFromMatrix(countData= countTable1,colData= colData2,design= ~Batch + Condition)
#This command line is to specify the reference group
dds2$Condition <- relevel(dds2$Condition, ref = "CONTROL")

#####Calling differential expression
dds2<-DESeq(dds2)
res2<-results(dds2)
summary (res2)
resultados2 <- as.data.frame(res2)
res3= as.data.frame(res2)
#filtering to get only significant genes
DE_genes_pvalue2= subset.data.frame(res3, padj<0.05 & abs(log2FoldChange)>1)
#saving output to a csv file
write.csv(DE_genes_pvalue2, file = DEGs_output_file)
#volcano plot
png(volcanoCOVSCA_output_file, units = "in", width = 5, height = 5, res = 300)
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res3$log2FoldChange, -log10(res3$pvalue))
plot(res3$log2FoldChange, -log10(res3$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=1.5)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res3$log2FoldChange) > 1 & res3$padj < alpha 
text(res3$log2FoldChange[gn.selected],
     -log10(res3$padj)[gn.selected],
     lab=rownames(res3)[gn.selected ], cex=0.4)
dev.off()

#A better volcano plot with another R package
library(EnhancedVolcano)
png(volcanoplotCAVSM1_output_file, units = "in", width = 8, height = 8, res = 600)
EnhancedVolcano(res3,
                lab = rownames(res3),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)

dev.off()
#getting the normalized counts
dds2 <- estimateSizeFactors(dds2)
sizeFactors(dds2)
normalized_counts <- counts(dds2, normalized=TRUE)
all_Z=t(scale(t(normalized_counts)))
normalized_counts_file <- construct_output_filename("normalized_counts.png", label)
write.csv(normalized_counts, file = normalized_counts_file)
#getting the PCA plot
pca <- prcomp(t(all_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=colData2$Condition)
pca.data
PCA_output_file <- construct_output_filename("PCAplot.png", label)
png(PCA_output_file, units="in", width=8, height=8, res=600)
ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("My PCA Graph")
dev.off()



