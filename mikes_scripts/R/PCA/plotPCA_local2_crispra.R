### Adapted Naky's PCA script for CRISPR screen ###


library(RColorBrewer)
library(ggplot2)
library(matrixStats)
library(edgeR)


setwd('C:/kutluaylab/data/CRISPRa_screen/mageck/count/count_total/') # with slash in the end
expname <- 'CRISPRa_HeLa_sort' # prefix for naming the output png
output_dir <- "C:/kutluaylab/data/CRISPRa_screen/PCA/" # with slash, output directory


raw_count <- read.delim(file = "count_total.count.txt")
data_sort <- raw_count[, seq(15, 23)]
row.names(data_sort) <- raw_count$sgRNA


data <- DGEList(counts = data_sort)
data <- calcNormFactors(data)


## MDS PLOT ##
colors <- brewer.pal(length(colnames(data)), name = "Set3")
names(colors) <- colnames(data)

# mds <- plotMDS(data, col = colors, pch = 18, cex = 2)
# legend <- legend('topright', title = "Experiment", legend = colnames(data),
#                  col = colors, cex = .75, pch = 18, xpd = T, bty = "n", inset = c(0, -.162))

# png(filename = file.path(output_dir, paste0(expname, "_MDS_1.png")), width = 2000, height = 1400, res = 300)
# print(mds)
# dev.off()

pdf(file = file.path(output_dir, paste0(expname, "_MDS_1.pdf")))
plotMDS(data, col = colors, pch = 18, cex = 2)
legend('topright', title = "Experiment", legend = colnames(data),
       col = colors, cex = .75, pch = 18, xpd = T, bty = "n", inset = c(0, -.162))
dev.off()





cpm <- cpm(data)


n <- 200
select <- order(rowVars(as.matrix(cpm)), decreasing = TRUE)
topNdata <- cpm[select[1:n], ]


logcpm <- log10(topNdata + 1)
# logcpm <- log10(cpm + 1)    # <<< for no top n selection


data_t <- as.data.frame(t(logcpm))

pca <- prcomp(data_t, center = T, scale = F) # main PCA calculation step
vars <- pca$sdev ^ 2
var_perc <- sprintf("%.2f", round((vars / sum(vars)) * 100, 2))
df <- data.frame(pca$x)


# df$sample <- row.names(df)
sample_label <- c("E484D-1", "AD-1", "E484D-2", "AD-2", "WT-1", "WT-2", "AD-4", "PRE-1", "PRE-2")

df$sample <- sample_label
# df$sample <- factor(sample_label, levels = sample_label)   # <<< to fix sample label order in plot (default would be alphabetical)




## PLOTTING PCA ##

colPalette <- c("#AEC7E8", "#98DF8A", "#FF9896", "#C49C94", "#C5B0D5", "#FFBB78", "#D62728", "#FF7F0E", "#9467BD", "#8C564B", "#E377C2", "#1F77B4", "#2CA02C", "#F7B6D2")

x_max <- max(df$PC1)
x_min <- min(df$PC1)
y_max <- max(df$PC2)
y_min <- min(df$PC2)
y.range <- y_max - y_min
x.range <- x_max - x_min


p <- (ggplot(df, aes(x=PC1, y=PC2, color=sample)) # PC1 vs PC2
  + geom_point(size=4)

  + scale_x_continuous(limits=c(x_min - x.range/10, x_max + x.range/10))
  + scale_y_continuous(limits=c(y_min - y.range/10, y_max + y.range/10))
  + scale_color_manual(values = colPalette, name = "Samples")

  + xlab(paste0("PC1: ", var_perc[1], "% var_perc")) # to plot PC1 vs. PC2
  + ylab(paste0("PC2: ", var_perc[2], "% var_perc"))

  + ggtitle(label = "Principal Component Analysis")
  + theme_bw()
  + theme(legend.position = 'right',
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(size = 10),
          aspect.ratio = 1)
)


png(filename = file.path(output_dir, paste0(expname, "_PCA_200.png")), width = 2000, height = 1400, res = 300)
print(p)
dev.off()
