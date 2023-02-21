# Introduction

This repository describes step-by-step RNA sequencing (RNA-seq) analysis conducted in zebrafish (_Danio rerio_) by Kumar _et al._ (2023), where zebrafish were exposed to muscle injury as follows: (i) single injury (SI), (ii) multiple injuries (MI), and (iii) multiple injuries at a later time point (day 5; D5MI). The transcriptomes from these conditions with three replicates were compared with the zebrafish transcriptome without any injury, serving as the control.

# RNA-seq steps
## Step 1: Obtaining count matrices

We employed [```salmon```](https://salmon.readthedocs.io/en/latest/) (v 1.9.0) for obtaining the count matrices from the ```_.fastq_``` files using [GRCz11 (GCA_000002035.4)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000002035.6/) with parameters ```--gcBias``` and ```--validateMappings```. We also generated ```decoys.txt``` as described in ```salmon``` [vignette](https://salmon.readthedocs.io/en/latest/salmon.html?highlight=decoy#preparing-transcriptome-indices-mapping-based-mode). The ```salmon``` script for obtaining count matrices can be found [here](/scripts/quantifier.sh).

## Step 2: ```tximport``` for ```DESeq2```

After quantification, we used [```tximport```](https://bioconductor.org/packages/release/bioc/html/tximport.html) (v 1.26.1) to import counts generated at step 1. For genome-wide annotation, we used [```org.Dr.eg.db```](https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html) for every sample that can be found in [```samples.txt```](samples.txt). See the workflow [here](scripts/salmon2tximport.r).

## Step 3: ```DESeq2``` object and global response

Using the ```txi``` object that includes all count matrices, we generated ```dds``` object using [```DESeq2``` package](http://master.bioconductor.org/packages/release/bioc/html/DESeq2.html) (v 1.38.3). Following filtering out genes with low counts, we visualised PCA plot, top 1000 variable genes and a pool of candidate genes as in violin plots with statistics and a heat map. For the script, click [here](scripts/deseq2visualisations.r).

## Step 4: Differentially expressed genes

Using this [script](scripts/deseq2results.r), we observed differentially expressed genes (DEGs) in each condition compared with control:
1. SI vs Control
2. MI vs Control
3. D5MI vs Control

## Step 5: Gene Ontology analysis

For over-representation and gene set enrichment analysis (GSEA), we utilised [```clusterProfiler```](http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) (v 4.6.0), [```enrichplot```](http://bioconductor.org/packages/release/bioc/html/enrichplot.html) (v 1.18.3) and [```DOSE```](http://bioconductor.org/packages/release/bioc/html/DOSE.html) (v 3.24.2) packages for [Gene Ontology](scripts/gene_ontology.r) analyses, as well as [```ReactomePA```](http://bioconductor.org/packages/release/bioc/html/ReactomePA.html) (v 1.42.0) for [Reactome](scripts/reactome.r) analysis.

Step 6: KEGG pathway analysis

Similar to Gene Ontology analysis, we also looked at [```KEGG```](https://www.kegg.jp) pathways for several candidate pathways. The script can be found [here](scripts/kegg.r).

