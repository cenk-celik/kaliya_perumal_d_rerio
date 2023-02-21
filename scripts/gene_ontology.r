library(clusterProfiler)
library(DOSE)
library(org.Dr.eg.db)

# `condition` parameter can be adjusted for other comparisons
condition <- 'D5MI'
control <- 'Ctrl'
comparison <- paste0('condition_', condition, '_vs_', control)
threshold <- log2(2)

dir.create(paste0(mainDir, 'plots/', comparison, '/Gene_ontology'), recursive = F)

# read the data
dat <- paste0('./results/results.', comparison, '.csv')

# Biological Pathways
ont <- 'BP' # other ontologies # MF, CC or ALL

d <- read.csv(dat)
geneList <- d[, 3] #log2FC
names(geneList) <- as.character(d[, 10]) # ENTREZIDs

#filtering
geneList <- na.omit(geneList)
geneList <- geneList[unique(names(geneList))]
geneList <- sort(geneList, decreasing = T)

gene <- names(geneList)[abs(geneList) > 2]
gene <- na.omit(gene)
gene <- gene[!duplicated(gene)]

setwd(paste0('plots/', comparison, '/Gene_ontology'))

# Over-representation
ego <-
  enrichGO(
    gene = gene,
    universe = names(geneList),
    OrgDb = org.Dr.eg.db,
    ont = ont,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T,
    minGSSize = 10,
    maxGSSize = 500
  )

library(dplyr)

filter(ego, p.adjust < 0.05, qvalue < 0.05)
mutate(ego, geneRatio = parse_ratio(GeneRatio)) %>%
  arrange(desc(geneRatio))
select(ego,-geneID) %>% head(n = 30)
ego_for_plot <-
  mutate(ego, richFactor = Count / as.numeric(sub('/\\d+', '', BgRatio)))

library(ggplot2)

pdf(paste0('./enrichGO.', ont, '.', comparison, '.pdf'))
dotplot(ego_for_plot, 
        showCategory = 15, 
        font.size = 10,
        orderBy = 'x') + 
  xlab('rich factor') + ylab(NULL) +
  ggtitle('Enriched Gene Ontology')
dev.off()

# GSEA
library(enrichplot)

gse <-
  gseGO(
    geneList = geneList,
    OrgDb = org.Dr.eg.db,
    ont = ont,
    minGSSize = 7,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    eps = 0
  )

gse_for_plot <- arrange(gse, abs(NES)) %>%
  group_by(sign(NES), .add = T)

pdf(paste0('./gene_set_enrichment.', ont, '.', comparison, '.pdf'))
dotplot(gse_for_plot, 
        showCategory = 10, 
        split = '.sign',
        font.size = 7,
        orderBy = 'x') + 
  facet_grid(.~.sign) +
  ylab(NULL)
dev.off()

setwd(mainDir)
