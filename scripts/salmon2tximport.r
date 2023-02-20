mainDir <- paste0(system("echo $HOME/", intern = T), "d_rerio/") # mainDir <- getwd() if running on Windows

library(tximport)
samples <- read.table('samples.txt', header = T, sep = ',')
files <- file.path(paste0(mainDir,'quants'), paste0(samples$run, '_quant'), 'quant.sf')

names(files) <- samples$run
all(file.exists(files))

library(readr)
library(AnnotationDbi)
library(org.Dr.eg.db)

edb <- org.Dr.eg.db

k <- keys(edb, keytype = 'ENSEMBLTRANS')
tx2gene <- AnnotationDbi::select(edb, k, "ENSEMBL", "ENSEMBLTRANS"); head(tx2gene)

txi <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = T)
