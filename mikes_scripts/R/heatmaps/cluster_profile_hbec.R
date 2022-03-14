# Time course trend profiles of different clusters, run after tc_cpm_heatmap code.

library(dplyr)
library(ggplot2)
library(RColorBrewer)


## Setup ##

# path to the /time_course_heatmap/ folder, created in the tc_cpm_heatmap code.
outputpath <- "C:/kutluaylab/data/Vero/riboseq/time_course_heatmap"
diffex_path <- file.path(outputpath, "diffex")
prefix <- "Vero_ribo"
outputname <- paste0(prefix, "_cluster_profile_test")

# number of clusters in the time course heatmap
num_cluster <- 10



## Reading the gene and cluster membership information from sig_cpm file ##
sig_cpm <- read.delim(file.path(outputpath, "sig_zscore", paste0(prefix, "_TCHeatmap_sig_cpm_cluster")))

gene_cluster <- na.omit(data.frame(sig_cpm$cluster, row.names = rownames(sig_cpm)))
colnames(gene_cluster) <- "cluster"
print(paste0("gene count in gene_cluster: ", nrow(gene_cluster)))



## Taking the intersect between diffex_genes_full table genes and genes in gene_cluster ##
gene_intersect <- rownames(gene_cluster)

for (file in dir(path = diffex_path)) {
  diffex <- read.delim(file.path(diffex_path, file))
  gene_intersect <- intersect(gene_intersect, rownames(diffex))
}

print(paste0("gene count present in gene_cluster and all diffex tables: ", length(gene_intersect)))



## Querying the logFC information from diffex_genes_full tables and collecting the information in a table called logFC_cluster_table ##

logFC_cluster_table <- NULL

for (file in dir(path = diffex_path)) {
  diffex <- read.delim(file.path(diffex_path, file), row.names = NULL)
  diffex <- rename(diffex, genes = row.names)
  gene_logFC <- diffex[match(gene_intersect, diffex$genes), ][, c("genes", "logFC")] # extracts genes and logFC information in the order of gene_intersect
  print(na.omit(nrow(gene_logFC)))

  hpi <- strsplit(file, split = "_")[[1]][1] # needs diffex file name to be in "[timepoint]hpi_..." format
  print(hpi)
  time <- as.integer(substr(hpi, start = 1, stop = nchar(hpi) - 3)) # getting rid of the "hpi"
  gene_logFC[, "time"] <- rep(time, times = nrow(gene_logFC))

  gene_logFC[, "cluster"] <- gene_cluster[gene_intersect, ]

  logFC_cluster_table <- rbind(logFC_cluster_table, gene_logFC)
}



## Creating a mock table and adding it to logFC_cluster_table as time = 0 ##
mock_table <- data.frame(gene_intersect,
                         rep(0, times = length(gene_intersect)), # logFC
                         as.integer(rep.int(0, times = length(gene_intersect))), # time
                         gene_cluster[gene_intersect, ])
colnames(mock_table) <- c("genes", "logFC", "time", "cluster")

logFC_cluster_table <- rbind(mock_table, logFC_cluster_table)

dir.create(path = file.path(outputpath, "cluster_profile"))
write.csv(logFC_cluster_table,
          file = file.path(outputpath, "cluster_profile",
                           paste0(outputname, "_logFC_cluster_table.csv")))




## Calculating the mean for each time-cluster ##
representative.profiles <- logFC_cluster_table %>%
  group_by(time, cluster) %>%
  summarize(logFC = mean(logFC))



## Calculating how many genes are in each cluster ##
num.in.cluster <- data.frame(table(gene_cluster[gene_intersect, ]))
names(num.in.cluster) <- c("cluster", "n")

cluster.labels <- paste0("Cluster ", num.in.cluster$cluster, " (n=", num.in.cluster$n, ")")
# names(cluster.labels) <- seq(num_cluster)
names(cluster.labels) <- num.in.cluster$cluster



## Plotting ##

colPalette <- c("#AEC7E87F", "#98DF8A7F", "#FF98967F",
                "#C49C947F", "#C5B0D57F", "#FFBB787F",
                "turquoise", "pink", "maroon", "wheat3")


profile_plots <- (ggplot(logFC_cluster_table, aes(x = as.numeric(time),
                                                  y = logFC,
                                                  group = interaction(genes, cluster),
                                                  color = as.factor(cluster),
                                                  facet = cluster))

  + geom_line(size = 1, alpha = 0.5)

  + geom_line(size = 1, data = representative.profiles,
              aes(x = as.numeric(time), y = logFC, group = as.factor(cluster)),
              color = "black")

  + geom_point(size = 1.5, data = representative.profiles,
               aes(x = as.numeric(time), y = logFC, group = NULL),
               color = "black")


  + scale_y_continuous(name = "log2(fold-change over mock)")
  + scale_x_continuous(name = "Hours post-infection",
                       # breaks = c(24, 48, 72, 96), # HBEC only
                       breaks = c(2, 6, 12, 24), # Vero only
                       expand = c(0.05, 0))

  # + facet_grid(cluster ~ ., nrow = 3, scales="free")
               #labeller = labeller(cluster = cluster.labels, size = 6)

  + facet_wrap(~ cluster, nrow = 4, scales = "free",
               labeller = labeller(cluster = cluster.labels))

  # + scale_color_brewer(type = "qual", palette = "Set3")

  + scale_color_manual(values = colPalette)

  + theme_classic()

  + theme(legend.position = "none",
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          # axis.text = element_text(size=8),
          text = element_text(size = 15),
          strip.background = element_blank())
)


png(file.path(outputpath, "cluster_profile", paste0(outputname, ".png")), width = 1500, height = 2000, res = 300)
print(profile_plots)
dev.off()

