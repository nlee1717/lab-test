### Adapted Naky's PCA script for CRISPR screen using MAGeCK counts ###


library(RColorBrewer)
library(ggplot2)
library(matrixStats)
library(edgeR)
library(stringr)


setwd('C:/kutluaylab/data/CRISPRa_screen/ng_screen/mageck/count/count_total/') # with slash in the end
expname <- 'CRISPRa_ng_sort_e275' # prefix for naming the output png
output_dir <- "C:/kutluaylab/data/CRISPRa_screen/ng_screen/PCA/" # with slash, output directory


raw_count <- read.delim(file = "count_total.count.txt")
# data_selected <- raw_count[, 3:13]
# pattern for all sorts: "e27[0-9]_[a-z][0-9]"
# pattern for E278: "e278_[a-z][0-9]"
data_selected <- raw_count[, str_detect(string = colnames(raw_count), pattern = "e275")]
row.names(data_selected) <- raw_count$sgRNA


data <- DGEList(counts = data_selected)
data <- calcNormFactors(data) # calcNormFactors is good when there is no DE among most genes
write.csv(data$samples, file = file.path(output_dir, paste0(expname, "_norm_factors.csv")))


## MDS PLOT ##
colors <- brewer.pal(ncol(data), name = "Set3")
names(colors) <- colnames(data)

# mds <- plotMDS(data, col = colors, pch = 18, cex = 2)
# legend <- legend('topright', title = "Experiment", legend = colnames(data),
#                  col = colors, cex = .75, pch = 18, xpd = T, bty = "n", inset = c(0, -.162))

# png(filename = file.path(output_dir, paste0(expname, "_MDS_1.png")), width = 2000, height = 1400, res = 300)
# print(mds)
# dev.off()

# check width, mar, inset to ensure proper legend position
pdf(file = file.path(output_dir, paste0(expname, "_MDS.pdf")), width = 8)
par(mar = c(5, 4, 4, 10))
plotMDS(data, col = colors, pch = 18, cex = 2)
legend("topright", title = "Experiment", legend = colnames(data),
       col = colors, cex = 1, pch = 18, xpd = NA, bty = "n", inset = c(-0.35, 0))
dev.off()





cpm <- cpm(data)
write.csv(cpm, file = file.path(output_dir, paste0(expname, "_cpm.csv")))


# expname <- 'CRISPRa_ng_sort_e278_top500' # prefix for naming the output png
# n <- 500
# select <- order(rowVars(as.matrix(cpm)), decreasing = TRUE)
# topNdata <- cpm[select[1:n], ]

topNdata <- cpm   # <<< for no top n selection


logcpm <- log10(topNdata + 1)


data_t <- as.data.frame(t(logcpm))

pca <- prcomp(data_t, center = T, scale = F) # main PCA calculation step
vars <- pca$sdev ^ 2
var_perc <- sprintf("%.2f", round((vars / sum(vars)) * 100, 2))
df <- data.frame(pca$x)


# df$sample <- row.names(df)
# sample_label <- c("E484D-1", "AD-1", "E484D-2", "AD-2", "WT-1", "WT-2", "AD-4", "PRE-1", "PRE-2")
# sample_label <- c("24h mock", "24h WT", "24h E484D", "72h mock", "72h WT", "72h E484D", "24h mock", "24h WT", "24h E484D", "72h mock", "72h WT", "72h E484D")
# reps <- c(rep("rep1", 6), rep("rep2", 6))

# df$sample <- sample_label
# df$sample <- factor(sample_label, levels = sample_label)   # <<< to fix sample label order in plot (default would be alphabetical)

# df$sample <- factor(sample_label, levels = c("24h mock", "24h WT", "24h E484D", "72h mock", "72h WT", "72h E484D"))
# df$reps <- reps

# sample <- sapply(strsplit(row.names(df), split = '_'), function(x) x[1])
# df$sample <- factor(sample, levels = c("uninfected", "wt", "P7", "delta", "mu", "lambda"))
# df$reps <- paste0("rep", sapply(strsplit(row.names(df),split = '_'), function(x) x[2]))



## PLOTTING PCA ##

colPalette <- c("#AEC7E8", "#98DF8A", "#FF9896", "#C49C94", "#C5B0D5", "#FFBB78", "#D62728", "#FF7F0E", "#9467BD", "#8C564B", "#E377C2", "#1F77B4", "#2CA02C", "#F7B6D2")


# PC1 vs PC2
x_max <- max(df$PC1)
x_min <- min(df$PC1)
y_max <- max(df$PC2)
y_min <- min(df$PC2)
y.range <- y_max - y_min
x.range <- x_max - x_min


p <- (ggplot(df, aes(x = PC1, y = PC2, color = row.names(df))) # PC1 vs PC2
  + geom_point(size = 4)

  + scale_x_continuous(limits = c(x_min - x.range / 10, x_max + x.range / 10))
  + scale_y_continuous(limits = c(y_min - y.range / 10, y_max + y.range / 10))
  # + scale_shape_discrete(name = "Replicates")
  + scale_color_manual(values = colPalette, name = "Experiment")

  + xlab(paste0("PC1: ", var_perc[1], "% var_perc")) # to plot PC1 vs. PC2
  + ylab(paste0("PC2: ", var_perc[2], "% var_perc"))

  + ggtitle(label = "Principal Component Analysis")
  + theme_bw()
  + theme(legend.position = 'right',
          # legend.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(size = 10),
          aspect.ratio = 1)
)


png(filename = file.path(output_dir, paste0(expname, "_PCA_all.png")), width = 2000, height = 1400, res = 300)
print(p)
dev.off()



# PC2 vs PC3
x_max <- max(df$PC2)
x_min <- min(df$PC2)
y_max <- max(df$PC3)
y_min <- min(df$PC3)
y.range <- y_max - y_min
x.range <- x_max - x_min


p2 <- (ggplot(df, aes(x = PC2, y = PC3, color = row.names(df))) # PC2 vs PC3
  + geom_point(size = 4)

  + scale_x_continuous(limits = c(x_min - x.range / 10, x_max + x.range / 10))
  + scale_y_continuous(limits = c(y_min - y.range / 10, y_max + y.range / 10))
  # + scale_shape_discrete(name = "Replicates")
  + scale_color_manual(values = colPalette, name = "Experiment")

  + xlab(paste0("PC2: ", var_perc[2], "% var_perc")) # to plot PC1 vs. PC2
  + ylab(paste0("PC3: ", var_perc[3], "% var_perc"))

  + ggtitle(label = "Principal Component Analysis")
  + theme_bw()
  + theme(legend.position = 'right',
          # legend.title = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(size = 10),
          aspect.ratio = 1)
)


png(filename = file.path(output_dir, paste0(expname, "_PCA_all_2_vs_3.png")), width = 2000, height = 1400, res = 300)
print(p2)
dev.off()

