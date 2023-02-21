library(ggplot2)
library(AnnotationDbi)
library(org.Dr.eg.db)
library(DESeq2)

dds <- readRDS('RDS/dds.rds')
vsd <- readRDS('RDS/vsd.rds')

##----comparison----
# Here we show only single injury vs control; however, other comparisons can also be done by
# changing `condition` parameter to either one of MI or D5MI

condition <- 'SI'
control <- 'Ctrl'
comparison <- paste0('condition_', condition, '_vs_', control)

# log2 fold change threshold for significance testing
threshold <- log2(2)

dir.create(paste0(mainDir, 'plots/', comparison), recursive = F)

# run DEG analysis
res <-
  results(
    dds,
    name = comparison,
    lfcThreshold = threshold,
    alpha = 0.05
  )

summary(res)
sum(res$padj < 0.1, na.rm = T)

# create several different annotation columns in res object for later
res$symbol <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(res),
    column = 'SYMBOL',
    keytype = 'ENSEMBL',
    multiVals = 'first'
  )

res$name <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(res),
    column = 'GENENAME',
    keytype = 'ENSEMBL',
    multiVals = 'first'
  )

res$entrez <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(res),
    column = 'ENTREZID',
    keytype = 'ENSEMBL',
    multiVals = 'first'
  )

res$GO <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(res),
    column = 'GO',
    keytype = 'ENSEMBL',
    multiVals = 'first'
  )

# order results
res <- res[order(res$pvalue),]

# write in a csv file
dir.create(paste0(mainDir, 'results/'), recursive = F)
write.csv(as.data.frame(res), file = paste0('./results/results.', comparison, '.csv'))

# only significant genes
ressig <- subset(res, padj < 0.1)
write.csv(as.data.frame(ressig), file = paste0('./results/results.', comparison, '.sig.csv'))

# Overall DEGs:
rownames(dds) <-
  mapIds(
    org.Dr.eg.db,
    keys = row.names(dds),
    keytype = 'ENSEMBL',
    column = 'SYMBOL',
    multiVals = 'first'
  )

library(EnhancedVolcano)

up_genes <- rownames(res[res$log2FoldChange > 0 & res$pvalue < 10e-5, ])
down_genes <- rownames(res[res$log2FoldChange < 0 & res$pvalue < 10e-5, ])

keyvals.colour <- ifelse(
  rownames(res) %in% up_genes, 'blue',
  ifelse(rownames(res) %in% down_genes, 'red3', 'black')
)

keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'N.S.'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Upregulated'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'Downregulated'

pdf(paste0('./plots/', comparison, '/volcano.', comparison, '.pdf'), width = 7, height = 5)
EnhancedVolcano(
  res,
  lab = res$symbol,
  selectLab = '',
  x = 'log2FoldChange',
  title = comparison,
  y = 'padj',
  legendPosition = 'right',
  FCcutoff = threshold,
  legendLabels = NULL,
  drawConnectors = F,
  boxedLabels = F,
  subtitle = '',
  caption = '',
  labSize = 3,
  pointSize = 1,
  colAlpha = 0.5,
  colCustom = keyvals.colour
) + theme_minimal()
dev.off()

sessionInfo()
