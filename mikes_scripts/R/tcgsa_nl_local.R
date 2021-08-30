library(edgeR)
library(TcGSA)
library(GSA)
library(org.Hs.eg.db)
library(stringr)

setwd("C:/kutluaylab/data/Vero/Vero_riboseq_counts") # count file names need to be in experiment group_condition_timepoint format


## Gene Set
gmt <- GSA.read.gmt("C:/kutluaylab/data/annotation_files/hallmark.gmt")

gtf_path <- "C:/kutluaylab/data/annotation_files/GRCh37.87_notab_geneonly.gtf"


## Design Matrix
saveMock <- function(cond,time) {
  if (cond=='mock') return(paste(cond,time,sep='-'))
  else return(time)
}


samples <- filelist <- list.files(pattern = '*hpi')
expers <- sapply(strsplit(filelist, '_'), function(x) x[1])
cond <- sapply(strsplit(filelist, '_'), function(x) x[2])
times <- sapply(strsplit(filelist, '_'), function(x) x[3])
times <- unname(mapply(saveMock, cond, times))

for (i in seq_along(times)) {
  if (grepl('mock', times[i])) {times[i] <- 0} # mock is treated as 0h
  else {times[i] <- gsub(pattern = "[^[:digit:].]", replacement = "", times[i])} # replaces the "hpi" in timepoints with ""
}

times <- as.numeric(times)

design <- data.frame(samples, expers, cond, times)


## Gene Expression (counts)
DG <- readDGE(design[, 1], header=FALSE)

if (grepl(pattern = "ENSG", rownames(DG$counts)[1])) { # updated 5/28
  gtf <- read.table(gtf_path, sep = " ")
  gene_names <- as.character(gtf[, 6])
  gene_names <- substr(gene_names, start = 1, stop = nchar(gene_names) - 1) # getting rid of the ";" at the end
  print(gene_names[1:10])
  rownames(DG) <- gene_names
}

# symbols <- mapIds(org.Hs.eg.db, keys = rownames(DG$counts), keytype = "ENSEMBL", column="SYMBOL")
# rownames(DG) <- symbols
# symbols_to_keep <- na.omit(symbols)
# DG <- DG[symbols_to_keep,]

keep <- rowSums(cpm(DG) > 1) >= 0.5 * ncol(DG) # filtering by expression: cpm > 1 in at least half of the samples
print(paste0("before filtering: ", nrow(DG)))
DG <- DG[keep, ,keep.lib.sizes=FALSE]
print(paste0("after filtering: ", nrow(DG)))

DG <- calcNormFactors(DG)
cpm <- cpm(DG, log = T, prior.count = 1) # added prior.count = 1 5/28 (default is 2)
# saveRDS(DG, "DG.rds")
# write.table(cpm, "cpm", quote=FALSE, row.names=TRUE, col.names=TRUE)


## TcGSA Model Fitting
tcgsa_res <- TcGSA::TcGSA.LR(expr = cpm,
                             gmt = gmt,
                             design = design,
                             subject_name = "expers",
                             time_name = "times")

summary(tcgsa_res)
head(TcGSA::signifLRT.TcGSA(tcgsa_res)$mixedLRTadjRes)


## Clustering
# clust <- TcGSA::clustTrend(tcgs = tcgsa_res, 
#                            expr =  tcgsa_res$Estimations,
#                            Subject_ID = design$expers,
#                            TimePoint = design$times)
clust <- TcGSA::clustTrend(tcgs = tcgsa_res, 
                           expr = cpm, # why use cpm here, instead of tcgsa_res$Estimations?
                           Subject_ID = design$expers,
                           TimePoint = design$times,
                           myproc = "HOLM", only.signif = TRUE) # why use "HOLM" method here?

print(clust)


## Heatmap
png('TcGSA_heatmap.png', height=2000, width=1000)
par(mar=c(5.1, 4.1, 4.1, 3.1))
plot(x = tcgsa_res, expr = cpm, clust_trends = clust,
     Subject_ID = design$expers, TimePoint = design$times,
     margins = c(10, 30), heatmap.width = 0.2, heatmap.height = 7,
     dendrogram.size = 0.05, heatKey.size = 0.5, cex.label.row = 4.1,
     cex.label.col = 4, myproc = "HOLM")
dev.off()
