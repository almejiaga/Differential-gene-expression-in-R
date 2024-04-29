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

![image](https://github.com/almejiaga/Differential-gene-expression-in-R/assets/124840761/669dee2f-c226-46f0-b67e-f793cf1d3b2a)

### 2. Normalized counts
"","64","66","71","73","78","80"
"ENSG00000268903.1",26.4855897629347,34.7521870914913,45.6143084145227,34.0329617965659,24.7432308631646,27.2958997615103
"ENSG00000225630.1",1797.96407429153,1482.01262370811,1505.27217767925,1426.38148706195,972.599305467471,1009.94829117588
"ENSG00000237973.1",479.796645320855,573.971606156243,377.016222609831,500.48473230244,275.03052767133,243.71339072777
"ENSG00000248527.1",4501.53158316955,5060.36685583844,3438.760311903,4462.32187320855,1734.88114859804,1819.07674839208
"ENSG00000198744.5",94.7369172289586,103.1355229812,162.908244337581,97.0940380666733,36.1631835692406,37.0444353906211
"ENSG00000228327.4",24.4482367042474,39.2363402645869,33.5125531208738,60.0581678762928,42.8248226477849,33.1450211389768
"ENSG00000228794.11",344.312666918151,464.109853415399,432.870477811287,513.497335342303,510.091220871394,459.15602813112
"ENSG00000230699.2",23.4295601749037,25.7838807453,58.6469679615292,64.0620457347123,51.3897871773419,39.9689960793544
"ENSG00000223764.2",24.4482367042474,45.9625700242304,42.8215956544499,82.0794960976001,50.4381244518356,56.5415066488427
"ENSG00000188976.11",2324.61983996219,2134.45691039353,1982.8260596517,2132.06495960839,2070.81809070178,1851.24691596814
"ENSG00000187961.15",326.995165919309,440.568049256647,397.496116183698,548.531266603474,500.57459361633,404.564228608099
"ENSG00000187583.11",265.87457415869,304.922415770504,289.511222794216,363.351915651571,318.807013044621,232.015147972837

### 3. Basic Volcano plot with enhanced volcano plot R package
![OBECO_enhancedvolcano](https://github.com/almejiaga/Differential-gene-expression-in-R/assets/124840761/bdeb8a79-8f78-4d7e-935a-d550a08fe1ef)

### 4. PCA plot with PC1 and PC2 colored by treatment

![OBECO_PCAplot](https://github.com/almejiaga/Differential-gene-expression-in-R/assets/124840761/58c4ce2a-6d55-4bbf-99da-267455f6bca5)

## Citation

If you use any of this code for your analysis, please cite:

Mejia-Garcia A, Fernandez GJ, Echeverri LF, Balcazar N, Acin S. RNA-seq analysis reveals modulation of inflammatory pathways by an enriched-triterpene natural extract in mouse and human macrophage cell lines. Heliyon. 2024 Jan 10;10(2):e24382. doi: 10.1016/j.heliyon.2024.e24382. PMID: 38293365; PMCID: PMC10826738.



##
