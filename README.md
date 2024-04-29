# Differential-gene-expression-in-R
Differential gene expression analysis in R using Deseq2 to normalize and scale data, generate PCA plots, Perform K means clustering, plotting Heatmaps, computing correlations, saving normalized matrices, perform annotation and GSEA.

## Input for the analysis

You will need two files:

A) Gene expression matrix. In this case, this input was obtained from FeatureCounts, and the first row corresponds to rownames and each column is an individual

```
			TCGA1		TCGA2		TCGA3		TCGA4	
ENSG00000000003.13	11.04916787	9.665335917	10.01262454	11.14784089	
ENSG00000000460.15	7.247927513	6.426264755	7.930737338	7.247927513	
ENSG00000000971.14	6.918863237	6.50779464	7.499845887	7.977279923	
ENSG00000001036.12	10.666224	9.826548487	12.07881795	11.07347215	
ENSG00000001084.9	9.30833903	9.157346935	9.812177306	11.11634396	
ENSG00000001167.13	11.26795708	10.67683861	10.71338651	11.86534664	
ENSG00000001460.16	9.853309555	7.794415866	8.550746785	9.519636253	
ENSG00000001461.15	10.69348696	8.864186145	9.618385502	9.958552715	
ENSG00000001497.15	11.38046107	10.97441459	11.62981194	10.97226185

```
B) A factors table where you store the treatment for each individual and other covariables that you would like to include in the analysis
```
sample	Condition	batch
1	CONTROL		exp1
2	CD		exp1
3	OBE		exp1
4	M1		exp1
5	UA		exp1
6	OA		exp1
7	UAL		exp1
8	CONTROL		exp2
```

## Running the scriopt
Rscript DEGs.R --count_table macrofagoshumanos_counts.txt --factors_table Factors.txt --condition1 "CONTROL" --condition2 "OBE" --label "OBECO"

Additional arguments:

--condition1: the first condition (its going to be the reference for the analysis

--condition2: the alternative condition

--label: label to include in all the outputs for the analysis

## Outputs 

### 1. DEGs csv file

This is the standard Differential gene expression output from Deseq2 but as a csv file

![image](https://github.com/almejiaga/Differential-gene-expression-in-R/assets/124840761/669dee2f-c226-46f0-b67e-f793cf1d3b2a)

### 2. Normalized counts

the normalized counts by depth and library composition as it is done by Deseq2, you can notice that now we have floats in the counts, indicating that they have been normalized


"","64","66"
"1",26.2,34.7
"2",1797.9,1482.01
"3",479.5,573.9
"4",4501.5,5060.36
"5",94.7,103.13


### 3. Basic Volcano plot with enhanced volcano plot R package

Includes color for significant genes, and their names (in this case Ensembl ID)

![OBECO_enhancedvolcano](https://github.com/almejiaga/Differential-gene-expression-in-R/assets/124840761/bdeb8a79-8f78-4d7e-935a-d550a08fe1ef)

### 4. PCA plot with PC1 and PC2 colored by treatment

PC1 vs PC2 with a different color by treatment

![OBECO_PCAplot](https://github.com/almejiaga/Differential-gene-expression-in-R/assets/124840761/58c4ce2a-6d55-4bbf-99da-267455f6bca5)

## Citation

If you use any of this code for your analysis, please cite:

Mejia-Garcia A, Fernandez GJ, Echeverri LF, Balcazar N, Acin S. RNA-seq analysis reveals modulation of inflammatory pathways by an enriched-triterpene natural extract in mouse and human macrophage cell lines. Heliyon. 2024 Jan 10;10(2):e24382. doi: 10.1016/j.heliyon.2024.e24382. PMID: 38293365; PMCID: PMC10826738.



##
