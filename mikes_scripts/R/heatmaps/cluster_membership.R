path <- "C:/kutluaylab/data/Vero/rnaseq/time_course_heatmap"
name <- "Vero_rseq"

sig_cpm <- read.delim(file.path(path, "sig_zscore", paste0(name, "_TCHeatmap_sig_cpm_cluster")))

gene_cluster <- na.omit(data.frame(row.names(sig_cpm), sig_cpm$cluster))
colnames(gene_cluster) <- c("gene_name", "cluster")

gene_cluster_sorted <- gene_cluster[order(gene_cluster$cluster), ]
row.names(gene_cluster_sorted) <- seq(nrow(gene_cluster_sorted))

write.csv(gene_cluster_sorted, file.path(path, paste0(name, "_cluster_membership.csv")))
