# ---
# title: "PCA"
# author: "NaKyung Lee"
# date: "3/8/2021"
# output: html_document
# ---

library(RColorBrewer)
library(ggplot2)
#library(ggrepel)
library(matrixStats)
library(edgeR)


# # Set directory to the folder with counts files
# # filename format should be: experiment_condition_timepoint
# # condition = infected OR mock

setwd('C:/kutluaylab/data/CRISPRa_screen/mageck/count/count_total/') # with slash in the end
expname <- 'CRISPRa_HeLa_sort' # prefix for naming the output png
output_dir <- "C:/kutluaylab/data/CRISPRa_screen/PCA/" # with slash, output directory

### LOADING DATA
filelist <- list.files(pattern = "*hpi") ## adjust the pattern accordingly
print(filelist)
data <- read.table(filelist[1], head = FALSE, sep = '\t')
row.names(data) <- data[,1]
data <- data[, -c(1, 2)]

for (file in filelist) {
  table <- read.table(file, head = FALSE, sep = '\t', col.names=c('id', file))
  data <- cbind(data, table[, 2])
}

colnames(data) <- filelist


### SETTING UP DESIGN MATRIX
conv_to_rep <- function(explist) { # converts the experiment names to "repn"
  exps <- unique(explist)
  #print(explist)
  for (i in seq_along(exps)) {
    explist <- gsub(exps[i], paste0('rep', toString(i)), explist)
  }
  #print(explist)
  return(explist)
}

saveMock <- function(cond,time) { # saves the time points of the mock samples
  if (cond == 'mock') return(paste(cond, time, sep = '-')) # saves as "mock-nhpi"
  #if (cond == 'mock') return(cond) # saves as "mock"
  else return(time)
}

experiments <- sapply(strsplit(filelist,'_'), function(x) x[1])
experiments <- conv_to_rep(experiments)
print(experiments)
condition <- sapply(strsplit(filelist,'_'), function(x) x[2])
print(condition)
timepoints <- sapply(strsplit(filelist,'_'), function(x) x[3])
timepoints <- unname(mapply(saveMock, condition, timepoints))
print(timepoints)

design <- data.frame(experiments, condition, timepoints)
rownames(design) <- filelist


### for CRISPRa screen PCA only ###
raw_count <- read.delim(file = "count_total.count.txt")
data_sort <- raw_count[, seq(15, 23)]
row.names(data_sort) <- raw_count$sgRNA
###################################


### PREPARING FOR PCA PLOT

## Data Normalization
data <- DGEList(counts = data, group = timepoints)
#data <- DGEList(counts = data, group = condition)
print(paste0("before filtering: ", nrow(data$counts)))
#print(data$samples$lib.size)

keep <- filterByExpr(data)
data <- data[keep, , keep.lib.sizes = F]
print(paste0("after filtering: ", nrow(data$counts)))

data <- calcNormFactors(data)
#print(data$samples$norm.factors)

data <- cpm(data)


## selecting genes with cpm > 1 across all samples

topNdata <- data[rowSums(data > 1) == ncol(data), ]
print(nrow(topNdata))
#print(topNdata[1:20, ])

## selecting top n genes with the most variance (needs at least n genes that pass the previous filters)

n <- 1000
select <- order(rowVars(as.matrix(topNdata)), decreasing = TRUE)
topNdata <- topNdata[select[1:n], ]


## transforming the normalized counts into log10 scale

topNdata <- log10(topNdata)
#topNdata <- log10(topNdata + 1)


## finding var_perc for PC1 and PC2 (or PC3 as well)

data_t <- as.data.frame(t(topNdata))
#data_t <- as.data.frame(t(data))

pca <- prcomp(data_t, center = T, scale = F) # main PCA calculation step
vars <- pca$sdev ^ 2
var_perc <- sprintf("%.2f", round((vars / sum(vars)) * 100, 2))
#print(var_perc)

df <- cbind(data.frame(pca$x), design)

## designing labels for the plot (ggrepel) (not currently used)
df$Label <- df$timepoints
#df$Label[which(design$experiments != 'rep1')] <- ""


### PLOTTING PCA

# Custom color palette
colPalette <- c("#AEC7E8", "#98DF8A", "#FF9896", "#C49C94", "#C5B0D5", "#FFBB78", "#D62728", "#FF7F0E", "#9467BD", "#8C564B", "#E377C2", "#1F77B4", "#2CA02C", "#F7B6D2")

## setting the x-axis and y-axis range
x_max <- max(df$PC1)
x_min <- min(df$PC1)
y_max <- max(df$PC2)
y_min <- min(df$PC2)
y.range <- y_max - y_min
x.range <- x_max - x_min

  df$timepoints <- factor(df$timepoints, levels = c("mock-2hpi", "mock-6hpi", "mock-12hpi", "mock-24hpi", "mock-hpi")) # to specify the ordering of timepoint legends for a particular plot

  p <- (ggplot(df, aes(x=PC1, y=PC2, shape=as.factor(experiments), color=timepoints, label=Label)) # PC1 vs PC2
        + geom_point(size=4)
        #+ geom_text_repel(size=4, box.padding=0.5)
        
        + scale_x_continuous(limits=c(x_min - x.range/10, x_max + x.range/10))
        + scale_y_continuous(limits=c(y_min - y.range/10, y_max + y.range/10))
        + scale_shape_discrete(name = "Experiment")
        + scale_color_manual(values = colPalette, name = "Timepoint")

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



## DIMNESION CHANGE (PC2 vs. PC3)
x_max <- max(df$PC2)
x_min <- min(df$PC2)
y_max <- max(df$PC3)
y_min <- min(df$PC3)
y.range <- y_max - y_min
x.range <- x_max - x_min

df$timepoints <- factor(df$timepoints, levels = c("mock-2hpi", "mock-6hpi", "mock-12hpi", "mock-24hpi", "mock-hpi")) # to specify the ordering of timepoint legends for a particular plot

p2 <- (ggplot(df, aes(x=PC2, y=PC3, shape=as.factor(experiments), color=timepoints, label=Label)) # PC2 vs PC3
        + geom_point(size=4)
        #+ geom_text_repel(size=4, box.padding=0.5)
        
        + scale_x_continuous(limits=c(x_min - x.range/10, x_max + x.range/10))
        + scale_y_continuous(limits=c(y_min - y.range/10, y_max + y.range/10))
        + scale_shape_discrete(name = "Experiment")
        + scale_color_manual(values = colPalette, name = "Timepoint")

        + xlab(paste0("PC2: ", var_perc[2], "% var_perc")) # to plot PC2 vs. PC3
        + ylab(paste0("PC3: ", var_perc[3], "% var_perc"))

        + ggtitle(label = "Principal Component Analysis")
        + theme_bw()
        + theme(legend.position = 'right',
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                panel.grid = element_blank(),
                plot.title = element_text(hjust=0.5),
                plot.caption = element_text(size = 10),
                aspect.ratio = 1)
  )


### SAVING THE PLOT
png(paste0(output_dir, paste(expname, "PCA_PC1_vs_PC2.png", sep = '_')),2000,1400, res = 300)
print(p)
dev.off()

png(paste0(output_dir, paste(expname, "PCA_PC2_vs_PC3.png", sep = '_')),2000,1400, res = 300)
print(p2)
dev.off()
