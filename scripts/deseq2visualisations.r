library(DESeq2)

# read sample information
coldata <- read.csv('samples.txt', header = T, sep = '\t')
colnames(txi$counts) <- coldata$names

# create dds object
dds <- DESeqDataSetFromTximport(txi,
                                colData = coldata,
                                design = ~condition)

# filter low counts
dds <- dds[rowSums(counts(dds)) > length(coldata$names), ]

# run DESeq
dds <- DESeq(dds)

# save for later
dir.create(paste0(mainDir, 'RDS/', recursive = F)
saveRDS(dds, file = 'RDS/dds.rds')
#dds <- readRDS('RDS/dds.rds')

# normalise
vsd <- vst(dds, blind = T)

# save for later
saveRDS(vsd, file = 'RDS/vsd.rds')
#vsd <- readRDS('RDS/vsd.rds')

# Principal component analysis
dir.create(paste0(mainDir, 'plots/', recursive = F)
pdf('./plots/pca.pdf')
plotPCA(vsd, intgroup = 'condition')
dev.off()

library(AnnotationDbi)
library(org.Dr.eg.db)

# replace ENSEBML gene names with SYMBOLs in vsd object
rownames(vsd) <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(vsd),
    keytype = 'ENSEMBL',
    column = 'SYMBOL',
    multiVals = 'first'
  )

# plot top 1000 variable genes heat map
topVarGenes <- head(order(-rowVars(assay(vsd))), 1000)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[ , 'condition'])
colnames(df) <- 'Treatment'
rownames(df) <- colnames(vsd)

library(pheatmap)

pdf('./plots/heatmap_top1000.pdf', width = 5, height = 7.5)
pheatmap(
  mat,
  cluster_rows = T,
  show_rownames = F,
  cluster_cols = T,
  annotation_col = df,
  angle_col = 45,
  color = viridis::inferno(100),
  cutree_rows = T,
  border_color = NA,
  scale = 'row',
  #cellwidth = 15,
  #cellheight = 10
)
dev.off()

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf('./plots/heatmap_sampleDistances.pdf')
pheatmap(sampleDistMatrix, cluster_rows = F, cluster_cols = F, angle_col = 45,
         annotation_col = df, color = viridis::inferno(1000), border_color = NA,
         cellwidth = 20, cellheight = 20)
dev.off()

##----Genes of interest----
library(ggpubr)
library(ggplot2)

# replace ENSEMBL gene names with SYMBOLs in dds object
rownames(dds) <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(dds),
    keytype = 'ENSEMBL',
    column = 'SYMBOL',
    multiVals = 'first'
  )

# list of genes to plot
selectedGenes <- c('alp3','col1a1a', #osteoblast diff-early
                   'runx2a','runx2b', #master switch
                   'bglap','pprc1','spp1','phex','sparc', #differentiation
                   'junba','smad1','smad5','ets1','cebpb','cebpg','dlx5a','her6','men1', #TF partners of runx2
                   'dlx3b','lef1', 'msx2a','msx2b','pparg','smad3a','hey1','stat1a', 'stat1b', #runx2 antogonist
                   'bmp2a','bmp3','bmp5','bmp6','bmp7b','hmox1a','fn1a',
                   'vegfaa','vegfab','vegfbb','vegfd',
                   'ctsk','runx3','lama1','scpp1','scpp5','scpp7','sall1a',
                   'sall1b','ocstamp','csf1ra','mfap4.12','mpeg1.3','mpeg1.2',
                   'shha','ptch2','smo','gli2a',
                   'dlx2a','dlx2b','dlx4a','dlx4b')
                   
# check if all gene names in dds object
selectedGenes %in% rownames(dds)

# the following function will generate violin plots for each gene in selectedGenes vector
lapply(
  selectedGenes,
  FUN = function(x) {
    expr <-
      plotCounts(
        dds,
        gene = x,
        intgroup = 'condition',
        returnData = T,
        normalized = T
      )
    pdf(paste0('./plots/violin_', x, '_expression.pdf'), height = 4, width = 6)
      ggplot(expr, aes(
        x = forcats::fct_relevel(condition, c('Ctrl', 'SI', 'MI', 'D5MI')),
        y = log2(count)
      )) +
      theme_minimal() + 
      labs(x = '', y = 'Count', title = x) +
      geom_violin(trim = F, aes(fill = condition)) +
      theme(legend.position = 'none') +
      geom_boxplot(width = 0.1, fill = 'transparent') +
      stat_compare_means(hide.ns = T,
        comparisons = list(c('SI', 'Ctrl'),
                           c('MI', 'Ctrl'),
                           c('D5MI', 'Ctrl')),
        label = 'p.format',
        ref.group = 'Ctrl',
        method = 't.test'
        )
    dev.off()
  }
)

# selected genes heat map
selectedGenesHeatmap <- c('alp3','col1a1a', #osteoblast diff-early
                   'runx2a','runx2b', #master switch
                   'bglap','pprc1','spp1','phex','sparc', #differentiation
                   'junba','smad1','smad5','ets1','cebpb','cebpg','dlx5a','men1', #TF partners of runx2
                   'dlx3b','lef1', 'msx2a','msx2b','smad3a','hey1','stat1a', 'stat1b', #runx2 antogonist
                   'bmp2a','bmp3','bmp5','bmp6','bmp7b','hmox1a','fn1a',
                   'vegfaa','vegfab','vegfbb','vegfd',
                   'ctsk','runx3','lama1','scpp1','scpp5','scpp7','sall1a',
                   'sall1b','ocstamp','csf1ra','mfap4.12','mpeg1.3','mpeg1.2',
                   'shha','ptch2','smo','gli2a',
                   'dlx2a','dlx2b','dlx4a','dlx4b'); selectedGenes %in% rownames(dds)


mat <- assay(vsd)[selectedGenesHeatmap, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[ , 'condition'])
colnames(df) <- 'Treatment'
rownames(df) <- colnames(vsd)

pdf('./plots/heatmap_selected.pdf', width = 5, height = 7.5)
pheatmap(
  mat,
  cluster_rows = T,
  show_rownames = T,
  cluster_cols = F,
  annotation_col = df,
  angle_col = 45,
  color = viridis::inferno(100),
  cutree_rows = T,
  border_color = NA,
  scale = 'row',
  cellwidth = 15,
  cellheight = 10
)
dev.off()
