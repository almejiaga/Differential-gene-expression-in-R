#Reading the gene expression matrix
countTable= read.table("macrofagoshumanos_counts.txt", header=T,row.names=1)
#removing unwanted columns from FeatureCounts output
countTable1 = countTable [,6:26]
countTable1 <- countTable1[,sort(colnames(countTable1))]
#organizing the names of the columns
vector4 <- seq(1:84)
colnames(countTable1) <- vector4[64:84]

#reading the factors table
colData= read.table ("Factors.txt", header=T)
#Filter to include only two groups for analysis
colData2 <- dplyr::filter(colData, Condition %in% c("CONTROL", "OBE"))
colData2$sample <- as.character(colData2$sample)
#filtering count table to exclude genes with less than 32 counts for any treatment
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
write.csv (DE_genes_pvalue2, file="DEGs_CO_vs_UAL.csv")
#volcano plot
png("volcanoCOVSCA.png", units="in", width=5, height=5, res=300)
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
png("volcanoplotCAVSM1.png", units="in", width=8, height=8, res=600)
EnhancedVolcano(res3,
                lab = rownames(res3),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0)

dev.off()
