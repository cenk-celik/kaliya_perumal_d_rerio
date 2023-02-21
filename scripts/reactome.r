library(ReactomePA)

# `condition` parameter can be changed for other comparisons
condition <- 'D5MI'
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

de <- names(geneList)[abs(geneList) > 2]
de <- na.omit(de)
de <- de[!duplicated(de)]

dir.create(paste0(mainDir, 'plots/', comparison, '/Reactome'), recursive = F)
setwd(paste0('plots/', comparison, '/Reactome'))

# Reactome analysis
library(ReactomePA)

x <- enrichPathway(
  gene = de,
  pvalueCutoff = 0.05,
  readable = T,
  organism = 'zebrafish'
)

filter(x, p.adjust < 0.05, qvalue < 0.2)

mutate(x, geneRatio = parse_ratio(GeneRatio)) %>%
  arrange(desc(geneRatio))

select(x, -geneID) %>% head

y <- mutate(x, richFactor = Count / as.numeric(sub('/\\d+', '', BgRatio)))

pdf(paste0('./enrichPathway.', '.', comparison, '.pdf'), width = 10, height = 5)
dotplot(y, 
        showCategory = 20, 
        font.size = 7,
        orderBy = 'x') + 
  xlab('rich factor') + 
  ylab(NULL) +
  ggtitle('Reactome Analysis')
dev.off()

#----Reactome GSEA----
y <-
  gsePathway(
    geneList,
    pvalueCutoff = 0.2,
    pAdjustMethod = 'BH',
    eps = 0,
    organism = 'zebrafish'
  )

t <- arrange(y, abs(NES)) %>%
  group_by(sign(NES), .add = T)

pdf(paste0('./gsePathway.', '.', comparison, '.pdf'), width = 10, height = 4)
dotplot(t, 
        showCategory = 10, 
        split = '.sign',
        font.size = 7,
        orderBy = 'x') + 
  facet_grid(.~.sign) +
  ylab(NULL)
dev.off()

setwd(mainDir)
