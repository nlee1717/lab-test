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

setwd('C:/kutluaylab/data/HBEC_162_163/rnaseq/counts/') # with slash in the end
expname <- 'hbec_rseq_test' # for naming the output png
output_dir <- "C:/kutluaylab/scripts/lab-test/mikes_scripts/R/PCA/test/" # with slash, output directory

### LOADING DATA
filelist <- list.files() ## adjust the pattern accordingly
print(filelist)
data <- read.table(filelist[1], head=FALSE, sep='\t')
row.names(data) <- data[,1]
data <- data[,-c(1,2)]

for (file in filelist) {
  table <- read.table(file, head=FALSE, sep='\t', col.names=c('id',file))
  data <- cbind(data, table[,2])
}

colnames(data) <- filelist


### SETTING UP DESIGN MATRIX
conv_to_rep <- function(explist) { # converts the experiment names to "repn"
  exps <- unique(explist)
  #print(explist)
  for (i in 1:length(exps)) {
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
condition <- sapply(strsplit(filelist,'_'), function(x) x[2])
print(condition)
timepoints <- sapply(strsplit(filelist,'_'), function(x) x[3])
timepoints <- unname(mapply(saveMock, condition, timepoints))
print(timepoints)

design <- data.frame(experiments, condition, timepoints)
rownames(design) <- filelist


### PREPARING FOR PCA PLOT

## Data Normalization
data <- DGEList(counts = data, group = timepoints)
#data <- DGEList(counts = data, group = condition)
print(nrow(data$counts))
#print(data$samples$lib.size)

keep <- filterByExpr(data)
data <- data[keep, , keep.lib.sizes=F]
print(nrow(data$counts))

data <- calcNormFactors(data)
#print(data$samples$norm.factors)

data <- cpm(data)
#data <- cpm(data, log = T, prior.count = 1)

#data <- log10(data + 1)


## selecting genes with cpm > 1 across all samples

topNdata <- data[rowSums(data > 1) == ncol(data), ]
print(nrow(topNdata))
#print(topNdata[1:20, ])

## selecting top n genes with the most variance

# n <- 1000
# select <- order(rowVars(as.matrix(topNdata)), decreasing = TRUE)
# topNdata <- topNdata[select[1:n],]


## transforming the normalized counts into log10 scale

topNdata <- log10(topNdata)
#topNdata <- log10(topNdata + 1)


## finding var_perc for PC1 and PC2 (or PC3 as well)

data_t <- as.data.frame(t(topNdata))
#data_t <- as.data.frame(t(data))

pca <- prcomp(data_t, center=T, scale=F) # main PCA calculation step
vars <- pca$sdev^2
var_perc <- sprintf("%.2f", round((vars/sum(vars)) * 100, 2))
#print(var_perc)

df <- cbind(data.frame(pca$x), design)

## setting the x-axis and y-axis range
#y.range <- max(df$PC2) - min(df$PC2)

min.x <- min(df$PC1) # to set PC1 as x axis
max.x <- max(df$PC1)

# min.x <- min(df$PC2) # to set PC2 as x axis (and PC3 as y axis)
# max.x <- max(df$PC2)

x.range <- max.x - min.x

## designing labels for the plot (ggrepel) (not currently used)
df$Label <- df$timepoints
df$Label[which(design$experiments != 'rep1')] <- ""

#df$condtime <- paste0(df$condition, ":", df$timepoints)


### PLOTTING PCA

  # Custom color palette
  colPalette <- c("#AEC7E8", "#98DF8A", "#FF9896", "#C49C94", "#C5B0D5", "#FFBB78", "#D62728", "#FF7F0E", "#9467BD", "#8C564B", "#E377C2", "#1F77B4", "#2CA02C", "#F7B6D2")

  #df$condtime <- factor(df$condtime, levels = c("CRISPR:6h", "CRISPR:24h", "CRISPR_IFN:6h", "CRISPR_IFN:24h", ))
  #df$timepoints <- factor(df$timepoints, levels = c("6h", "24h")) # to specify the ordering of timepoint legends for a particular plot
  df$timepoints <- factor(df$timepoints, levels = c("4hpi", "24hpi", "48hpi", "72hpi", "96hpi", "mock-4h", "mock-96h"))

  p <- (ggplot(df, aes(x=PC1, y=PC2, shape=as.factor(experiments), color=timepoints, label=Label)) # x=PC1, y=PC2 normally
        + geom_point(size=4)
        #+ geom_text_repel(size=4, box.padding=0.5)
        
        + scale_x_continuous(limits=c(min.x - x.range/10, max.x + x.range/10))
        + scale_shape_discrete(name = "Experiment")
        #+ scale_color_brewer(guide='none', type="qual", palette="Set2")
        + scale_color_manual(values = colPalette, name = "Timepoint")

        + xlab(paste0("PC1: ", var_perc[1], "% var_perc")) # to plot PC1 vs. PC2
        + ylab(paste0("PC2: ", var_perc[2], "% var_perc"))
        # + xlab(paste("PC2: ", var_perc[2], "% var_perc", sep="")) # to plot PC2 vs. PC3
        # + ylab(paste("PC3: ", var_perc[3], "% var_perc", sep=""))

        + ggtitle(label = "Principal Component Analysis")
        #+ labs(caption = "Note: Mock 24hpi was used for plotting.")
        + theme_bw()
        + theme(legend.position = 'right',
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                panel.grid = element_blank(),
                plot.title = element_text(hjust=0.5),
                #plot.subtitle = element_text(size=7, hjust=0.5, face="italic"),
                #plot.caption = element_text(size = 10),
                aspect.ratio = 1)
  )


### SAVING THE PLOT
#png(paste(expname,"PCA.png",sep='_'),350,350)
png(paste0(output_dir, paste(expname, "PCA.png", sep = '_')),2000,1400, res = 300)
p
dev.off()
