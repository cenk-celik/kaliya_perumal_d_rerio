# `condition` parameter can be replaced for other comparisons
condition <- 'SI'
control <- 'Ctrl'
comparison <- paste0('condition_', condition, '_vs_', control)
threshold <- log2(2)

dat <- paste0('./results/results.', comparison, '.csv')

d <- read.csv(dat)
geneList <- d[, 3]
names(geneList) <- as.character(d[, 10])
geneList <- na.omit(geneList)
geneList <- geneList[unique(names(geneList))]
geneList <- sort(geneList, decreasing = T)

gene <- names(geneList)[abs(geneList) > 0.5]
gene <- na.omit(gene)
gene <- gene[!duplicated(gene)]

dir.create(paste0(mainDir, 'plots/', comparison, '/kegg'), recursive = F)
setwd(paste0('plots/', comparison, '/kegg'))

dre <- search_kegg_organism('Danio rerio', by = 'scientific_name')

# KEGG pathway analysis
kk <- enrichKEGG(gene = gene,
                 organism = 'dre',
                 pvalueCutoff = 0.05)

library(pathview)

# only selected KEGG pathways to plot
path.ids <- c('dre04340', #Hedgehog
              'dre04310', #Wnt
              'dre04150', #mTOR
              'dre04330', #Notch
              'dre04350', #BMP
              'dre04621', #Nlrp3
              'dre04210', #Apoptosis
              'dre04020') #Ca

lapply(
  X = path.ids,
  FUN = function(x) {
    x <- pathview(
      gene.data = geneList,
      pathway.id = x,
      species = 'dre',
      limit = c(-5,5))
  }
)

# KEG GSEA for dotplot
kk2 <-
  gseKEGG(
    geneList = geneList,
    organism = 'dre',
    minGSSize = 10,
    pvalueCutoff = 0.05,
    eps = 0
  )

library(enrichplot)
library(ggplot2)

t <- arrange(kk2, abs(NES)) %>%
  group_by(sign(NES), .add = T)

pdf(paste0(getwd(),'/kegg_enrichment_dotplot.pdf'), width = 6,height = 10)
dotplot(t, 
        showCategory = 20, 
        split = '.sign',
        font.size = 7,
        orderBy = 'x') + 
  facet_grid(.~.sign) +
  ylab(NULL)
dev.off()

setwd(mainDir)
