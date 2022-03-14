### Time Course Heatmap ###


library(stringr)
library(edgeR)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


path <- "C:/kutluaylab/data/HBEC_162_163/rnaseq/time_course_heatmap/"
counts_path <- file.path(path, "counts")
diffex_path <- file.path(path, "diffex")
gtf_path <- "C:/kutluaylab/data/annotation_files/GRCh37.87_notab_geneonly.gtf"
outputname <- "HBEC_rseq_TCHeatmap_new"
numrep <- 4



## Reading Counts Data ##
filelist <- dir(path = counts_path, pattern = "16.+") # change accordingly
filelist <- str_sort(filelist, numeric = T)

data <- read.delim(file.path(counts_path, filelist[1]), header = F)
rownames(data) <- data[, 1]
data <- data[, -c(1, 2)]

for (file in filelist) {
  counts <- read.delim(file.path(counts_path, file), header = F)
  data <- cbind(data, counts[, 2])
}

colnames(data) <- filelist


saveMock <- function(cond, time) { # saves the time points of the mock samples
  #if (cond == 'mock') return(paste(cond, time, sep = '-')) # saves as "mock-nhpi"
  if (cond == 'mock') return(cond) # saves as "mock"
  else return(time)
}

experiment <- sapply(strsplit(filelist, '_'), function(x) x[1])
condition <- sapply(strsplit(filelist,'_'), function(x) x[2])
timepoints <- sapply(strsplit(filelist,'_'), function(x) x[3])
timepoints <- unname(mapply(saveMock, condition, timepoints))

design <- data.frame(experiment, condition, timepoints) # for checking
row.names(design) <- filelist



## DGE object, Filtering, CPM of a combined count table ##
dge <- DGEList(counts = data, samples = design, group = timepoints)
colnames(dge) <- paste0(experiment, "_", timepoints)
print(paste0("combined table before filtering: ", nrow(dge$counts)))


# Converting ensembl gene IDs to gene names using gtf file
if (grepl(pattern = "ENSG", rownames(dge$counts)[1])) {
  gtf <- read.table(gtf_path, sep = " ")
  gene_names <- as.character(gtf[, 6])
  gene_names <- substr(gene_names, start = 1, stop = nchar(gene_names) - 1) # getting rid of the ";" at the end
  rownames(dge) <- gene_names
}


# Filtering
keep <- rowSums(cpm(dge) > 1) >= numrep
dge <- dge[keep, ,keep.lib.sizes = F]
print(paste0("combined table after filtering: ", nrow(dge$counts)))

duplicate_genes <- row.names(dge)[duplicated(row.names(dge))]
dge <- dge[!(row.names(dge) %in% duplicate_genes), ,keep.lib.sizes = F]
print(paste0(length(duplicate_genes) + length(unique(duplicate_genes)),
             " duplicated genes removed."))

dge <- calcNormFactors(dge)
cpm <- cpm(dge, log = T, prior.count = 1)


write.table(cpm, file = file.path(path, "sig_zscore", paste0(outputname, "_cpm")), quote = F, sep = "\t")



## Selecting Significant DE Genes from edgeR Results and Taking Union ##
union_genes <- NULL

for (file in dir(path = diffex_path, pattern = "*full")) { # <<< change pattern
  diffex <- read.delim(file.path(diffex_path, file))
  sig <- rownames(diffex[abs(diffex[, "logFC"]) > 1 & diffex[, "FDR"] < 0.05, ]) # significant DE gene selection
  union_genes <- union(union_genes, sig)
}

print(paste0("union genes count: ", length(union_genes)))



## Significant CPM Table, Z-score ##
sig_cpm <- cpm[intersect(union_genes, row.names(cpm)), ]
print(paste0("union genes present in the filtered combined table: ", nrow(sig_cpm)))


# HBEC
col_order <- c("162-1_mock", "162-2_mock", "163-1_mock",
               "162-1_24hpi", "162-2_24hpi", "163-1_24hpi", "163-2_24hpi",
               "162-1_48hpi", "162-2_48hpi", "163-1_48hpi", "163-2_48hpi",
               "162-1_72hpi", "162-2_72hpi", "163-1_72hpi", "163-2_72hpi",
               "162-1_96hpi", "162-2_96hpi", "163-1_96hpi", "163-2_96hpi")

# Vero
# col_order <- c("176_mock", "rp3_mock", "rp9_mock",
#                "176_2hpi", "rp3_2hpi", "rp9_2hpi",
#                "176_6hpi", "rp3_6hpi", "rp9_6hpi",
#                "176_12hpi", "rp3_12hpi", "rp9_12hpi",
#                "176_24hpi", "rp3_24hpi", "rp9_24hpi")



sig_cpm <- sig_cpm[, col_order]
write.table(sig_cpm,
            file = file.path(path, "sig_zscore", paste0(outputname, "_sig_cpm")),
            quote = F,
            sep = "\t")

# for each gene as a unit, calculates the z-score of each sample
sig_zscore <- t(scale(t(as.matrix(sig_cpm)))) # z-score calculation
write.table(sig_zscore,
            file = file.path(path, "sig_zscore", paste0(outputname, "_sig_zscore")),
            quote = F,
            sep = "\t")




## Clustering ##
cluster <- ConsensusClusterPlus(t(sig_zscore),
                                maxK = 13,
                                reps = 100,
                                pItem = 0.8,
                                pFeature = 1,
                                # clusterAlg = "km",
                                # distance = "euclidean",
                                clusterAlg = "hc",
                                innerLinkage = "complete",
                                finalLinkage = "ward.D2",
                                distance = "pearson",
                                plot = "pdf",
                                title = file.path(path, "cluster", paste0(outputname, "_cluster"))
                                # seed = 1262118388.71284
                                )

save(cluster, file = file.path(path, "cluster", paste0(outputname, "_cluster.RData")))



## recreating heatmap starting from clustered RData and sig_zscore ##
path <- "C:/kutluaylab/data/Vero/rnaseq/time_course_heatmap/"
outputname <- "Vero_rseq_TCHeatmap"

load(file = file.path(path, "cluster", paste0(outputname, "_cluster.RData")))

sig_zscore <- read.delim(file = file.path(path, "sig_zscore", paste0(outputname, "_sig_zscore")))
sig_zscore_pruned <- sig_zscore



## Merging and Reordering Clusters ##
num.clusters <- 12
row.clusters.orig <- cluster[[num.clusters]][["consensusClass"]]
row.clusters <- row.clusters.orig
pos.1 <- which(row.clusters == 1)
pos.2 <- which(row.clusters == 2)
pos.3 <- which(row.clusters == 3)
pos.4 <- which(row.clusters == 4)
pos.5 <- which(row.clusters == 5)
pos.6 <- which(row.clusters == 6)
pos.7 <- which(row.clusters == 7)
pos.8 <- which(row.clusters == 8)
pos.9 <- which(row.clusters == 9)
pos.10 <- which(row.clusters == 10)
pos.11 <- which(row.clusters == 11)
pos.12 <- which(row.clusters == 12)


#HBEC_rseq version 1 (the one that was used)
# row.clusters[pos.1] <- 2
# row.clusters[pos.2] <- 2
# row.clusters[pos.3] <- 5
# row.clusters[pos.4] <- 99
# row.clusters[pos.5] <- 2
# row.clusters[pos.6] <- 3
# row.clusters[pos.7] <- 5
# row.clusters[pos.8] <- 4
# row.clusters[pos.9] <- 1
# row.clusters[pos.10] <- 2
# row.clusters[pos.11] <- 5
# row.clusters[pos.12] <- 6
#
# row.clusters <- row.clusters[-pos.4] # removing cluster 4
# sig_zscore_pruned <- sig_zscore[!(row.names(sig_zscore) %in% names(pos.4)), ] # removing cluster 4 genes from sig_zscore
# num.clusters <- 6

# HEBC rseq version 2
# row.clusters[pos.1] <- 2
# row.clusters[pos.2] <- 2
# row.clusters[pos.3] <- 5
# row.clusters[pos.4] <- 99
# row.clusters[pos.5] <- 2
# row.clusters[pos.6] <- 3
# row.clusters[pos.7] <- 5
# row.clusters[pos.8] <- 4
# row.clusters[pos.9] <- 1
# row.clusters[pos.10] <- 2
# row.clusters[pos.11] <- 5
# row.clusters[pos.12] <- 4
#
# row.clusters <- row.clusters[-pos.4] # removing cluster 4
# sig_zscore_pruned <- sig_zscore[!(row.names(sig_zscore) %in% names(pos.4)), ] # removing cluster 4 genes from sig_zscore
# num.clusters <- 5


#HBEC_ribo_version_1
# row.clusters[pos.1] <- 3
# row.clusters[pos.2] <- 2
# row.clusters[pos.3] <- 99
# row.clusters[pos.4] <- 5
# row.clusters[pos.5] <- 2
# row.clusters[pos.6] <- 4
# row.clusters[pos.7] <- 2
# row.clusters[pos.8] <- 2
# row.clusters[pos.9] <- 2
# row.clusters[pos.10] <- 6
# row.clusters[pos.11] <- 1
# row.clusters[pos.12] <- 6
#
# row.clusters <- row.clusters[-pos.3] # removing cluster 3
# sig_zscore_pruned <- sig_zscore[!(row.names(sig_zscore) %in% names(pos.3)), ] # removing cluster 3 genes from sig_zscore


#HBEC_ribo_version_2
# row.clusters[pos.1] <- 3
# row.clusters[pos.2] <- 2
# row.clusters[pos.3] <- 99
# row.clusters[pos.4] <- 5
# row.clusters[pos.5] <- 2
# row.clusters[pos.6] <- 4
# row.clusters[pos.7] <- 2
# row.clusters[pos.8] <- 2
# row.clusters[pos.9] <- 2
# row.clusters[pos.10] <- 5
# row.clusters[pos.11] <- 1
# row.clusters[pos.12] <- 5
#
# row.clusters <- row.clusters[-pos.3] # removing cluster 3
# sig_zscore_pruned <- sig_zscore[!(row.names(sig_zscore) %in% names(pos.3)), ] # removing cluster 3 genes from sig_zscore
# num.clusters <- 5


# Vero rnaseq
row.clusters[pos.1] <- 3
row.clusters[pos.2] <- 3
row.clusters[pos.3] <- 3
row.clusters[pos.4] <- 1
row.clusters[pos.5] <- 2
row.clusters[pos.6] <- 2
row.clusters[pos.7] <- 4
row.clusters[pos.8] <- 1
row.clusters[pos.9] <- 5
row.clusters[pos.10] <- 5
row.clusters[pos.11] <- 5
row.clusters[pos.12] <- 5

num.clusters <- 5

# Vero rnaseq test version
# row.clusters[pos.1] <- 3
# row.clusters[pos.2] <- 3
# row.clusters[pos.3] <- 3
# row.clusters[pos.4] <- 1
# row.clusters[pos.5] <- 2
# row.clusters[pos.6] <- 2
# row.clusters[pos.7] <- 4
# row.clusters[pos.8] <- 1
# row.clusters[pos.9] <- 5
# row.clusters[pos.10] <- 6
# row.clusters[pos.11] <- 7
# row.clusters[pos.12] <- 8
#
# num.clusters <- 8

# gene_cluster <- data.frame(names(row.clusters), as.integer(row.clusters))
# colnames(gene_cluster) <- c("gene_name", "cluster")
# row.names(gene_cluster) <- seq(nrow(gene_cluster))
# write.csv(gene_cluster, file.path(path, "Vero_rseq_cluster_membership_test.csv"))


# Vero riboseq
# row.clusters[pos.1] <- 5
# row.clusters[pos.2] <- 2
# row.clusters[pos.3] <- 8
# row.clusters[pos.4] <- 9
# row.clusters[pos.5] <- 3
# row.clusters[pos.6] <- 1
# row.clusters[pos.7] <- 3
# row.clusters[pos.8] <- 1
# row.clusters[pos.9] <- 10
# row.clusters[pos.10] <- 4
# row.clusters[pos.11] <- 7
# row.clusters[pos.12] <- 6
#
# num.clusters <- 10




sig_cpm <- read.delim(file = file.path(path, "sig_zscore", paste0(outputname, "_sig_cpm")))

outputname <- "HBEC_rseq_TCHeatmap_1"

write.table(sig_zscore_pruned,
            file = file.path(path, "sig_zscore", paste0(outputname, "_sig_zscore_pruned")),
            quote = F,
            sep = "\t")



## Saving cpm table of significant genes and their cluster membership info ##
sig_cpm <- data.frame(sig_cpm)
sig_cpm$cluster <- as.integer(row.clusters[match(rownames(sig_cpm), names(row.clusters))])
write.table(sig_cpm,
            file = file.path(path, "sig_zscore", paste0(outputname, "_sig_cpm_cluster")),
            quote = F,
            sep = "\t")


## Generating and saving cluster membership ##
gene_cluster <- na.omit(data.frame(row.names(sig_cpm), sig_cpm$cluster))
colnames(gene_cluster) <- c("gene_name", "cluster")

# generating cluster membership from row.clusters
gene_cluster <- data.frame(names(row.clusters), as.integer(row.clusters))
colnames(gene_cluster) <- c("gene_name", "cluster")

gene_cluster_sorted <- gene_cluster[order(gene_cluster$cluster), ]
row.names(gene_cluster_sorted) <- seq(nrow(gene_cluster_sorted))

write.csv(gene_cluster_sorted, file.path(path, paste0(outputname, "_cluster_membership.csv")))





### Plotting Heatmap ###

path <- "C:/kutluaylab/data/Vero/riboseq/time_course_heatmap/"
exp <- "Vero_ribo"
# outputname <- "HBEC_rseq_TCHeatmap_plotmod_test"


## Taking sig_zscore_pruned and cluster_membership as input ##

sig_zscore_pruned <- read.delim(file = file.path(path, "sig_zscore", paste0(exp, "_TCHeatmap_sig_zscore"))) # Vero
sig_zscore_pruned <- read.delim(file = file.path(path, "sig_zscore", paste0(exp, "_TCHeatmap_sig_zscore_pruned"))) # HBEC
# getting rid of the the "X" in front of the colnames (HBEC)
colnames(sig_zscore_pruned) <- substr(colnames(sig_zscore_pruned), start = 2, stop = nchar(colnames(sig_zscore_pruned)))

gene_cluster_sorted <- read.csv(file = file.path(path, paste0(exp, "_cluster_membership.csv")), row.names = 1)

# re-creating variable row.clusters
row.clusters <- gene_cluster_sorted$cluster
names(row.clusters) <- gene_cluster_sorted$gene_name

stopifnot(nrow(sig_zscore_pruned) == length(row.clusters) &
            setequal(row.names(sig_zscore_pruned), names(row.clusters)))
row.clusters <- row.clusters[row.names(sig_zscore_pruned)]


# Vero only, get rid of "X" in front of "176" #
colnames(sig_zscore_pruned) <- c("176_mock", "rp3_mock", "rp9_mock",
                                 "176_2hpi", "rp3_2hpi", "rp9_2hpi",
                                 "176_6hpi", "rp3_6hpi", "rp9_6hpi",
                                 "176_12hpi", "rp3_12hpi", "rp9_12hpi",
                                 "176_24hpi", "rp3_24hpi", "rp9_24hpi")



## ISGs to be labeled in the heatmap ##

# HBEC
isg <- c("IFI6", "XAF1", "IFITM1", "IFI44L", "AIM2", "OAS1",
         "IFI44", "CCL20", "IFIT1", "OAS2", "OAS3", "IFIH1",
         "IFIT3", "CMPK2", "OASL", "IDO1", "ISG15", "IFNL1",
         "MX2", "IFIT5", "IL6", "CXCL11", "CXCL12", "CXCL3",
         "CXCL6", "CXCL9", "CXCL2", "SELL", "CARD17", "NOS2", "TNFAIP6")
isg_tc <- isg[isg %in% names(row.clusters)]
isg_2 <- isg[isg %in% names(row.clusters[row.clusters == 2])] # the immune response cluster
gene_label_pos <- which(row.names(sig_zscore_pruned) %in% isg_2)
gene_label <- row.names(sig_zscore_pruned)[gene_label_pos]

# Vero
isg <- c("NR4A3", "OASL", "FOS", "CXCL8", "CXCL11", "NFKBIA",
         "CXCL3", "IFI44", "ISG15", "EGR1", "EGR2", "IFIT1",
         "IFIT2", "IFIT3", "CXCR4", "IL6", "IRF1", "IL11")
isg_tc <- isg[isg %in% names(row.clusters)]
# Vero
isg_1 <- isg[isg %in% names(row.clusters[row.clusters == 1])]
isg_2 <- isg[isg %in% names(row.clusters[row.clusters == 2])]
isg_12 <- isg[isg %in% names(row.clusters[row.clusters == 1 | row.clusters == 2])] # both cluster 1 and 3 have ISGs
gene_label_pos <- which(row.names(sig_zscore_pruned) %in% isg_12)
gene_label <- row.names(sig_zscore_pruned)[gene_label_pos]



# HBEC_162_163
split_order <- c("a", "a", "b", "a", "a", "b", "b", "a", "a", "b", "b", "a", "a", "b", "b", "a", "a", "b", "b")

# Vero_176_RP3_RP9
split_order <- rep(c("a", "b", "c"), 5)


colPalette <- c("#AEC7E87F", "#98DF8A7F", "#FF98967F",
                "#C49C947F", "#C5B0D57F", "#FFBB787F",
                "paleturquoise", "seagreen1", "orchid1", "wheat3") # custom color palette


print(max(row.clusters))
block.labels <- c("Immune Res.       \nClusters       ", "", "", "", "", "", "", "", "", "")
cluster.rowAnnot <- rowAnnotation(block = anno_block(gp = gpar(fill = colPalette,
                                                               col = NA),
                                                     labels = block.labels,
                                                     labels_gp = gpar(fontsize = 8)),
                                  width = unit(4, "mm"))
                                  #space = unit(2, "mm"))


gene.rowAnnot <- rowAnnotation(gene.name = anno_mark(at = gene_label_pos,
                                                     labels = gene_label,
                                                     labels_gp = gpar(fontsize = 7),
                                                     link_gp = gpar(lwd = 0.5),
                                                     link_width = unit(4, "mm"),
                                                     padding = unit(0.5, "mm")))


col_fun <- colorRamp2(breaks = seq(-3, 3, length = 256), # <<< can change data range for colors
                      colors = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)))


heatmap <- Heatmap(as.matrix(sig_zscore_pruned),
                   # labels
                   column_title = paste0("Differentially Expressed Genes\nn = ", nrow(sig_zscore_pruned)),
                   column_title_gp = gpar(fontsize = 9),
                   # column_order = col_order,
                   show_row_names = F,
                   column_names_gp = gpar(fontsize = 7),
                   name = "z-score",
                   col = col_fun,

                   # legends
                   show_heatmap_legend = F,
                   heatmap_legend_param = list(color_bar = "continuous",
                                               title_gp = gpar(fontsize = 6),
                                               labels_gp = gpar(fontsize = 6),
                                               grid_width = unit(2, units = "mm")), # not currently used
                   # clustering
                   cluster_columns = F,
                   clustering_distance_rows = "pearson",
                   cluster_row_slices = F,
                   show_row_dend = F,
                   split = row.clusters,
                   # splitting
                   column_split = factor(split_order),
                   cluster_column_slices = F,
                   column_gap = unit(2, units = "mm"),
                   # labels
                   # row_title_rot = 0,
                   # row_title_gp = gpar(fontsize = 7), # is this used?
                   row_title = NULL,
                   # annotation
                   left_annotation = cluster.rowAnnot,
                   right_annotation = gene.rowAnnot,
                   # size
                   width = ncol(sig_zscore_pruned) * 0.3,
                   # disable raster graphics
                   use_raster = F
)


# custom color legend
color_legend <- Legend(col_fun = col_fun,
                       title = "z-score",
                       title_gp = gpar(fontsize = 7),
                       title_position = "topcenter",
                       labels_gp = gpar(fontsize = 7),
                       grid_width = unit(2, units = "mm"))


outputname <- paste0(exp, "_TCHeatmap_test")

png(file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1000, res = 300)
draw(heatmap)
draw(color_legend, x = unit(0.9, "npc"), y = unit(0.3, "npc"))
dev.off()

