args <- commandArgs()

directory <- args[6]
files <- unlist(strsplit(args[7],","))
donor_list <- unlist(strsplit(args[9],","))
conditions <- unlist(strsplit(args[8],","))
n <- args[10]
prefix <- args[11]
gtf_path <- args[12]

print(files)
print(donor_list)
print(conditions)

library(edgeR)
library(statmod)
library(RColorBrewer)
library(gplots)
# library(org.Hs.eg.db)
library(ComplexHeatmap)
# library(DESeq2)
library(ggplot2)



#files <- c(filelist); #filelist is a list of comma separated file names
setwd(directory)
labels <- paste(donor_list, conditions, sep="_")
print(labels)
DG <- readDGE(files, header=FALSE, labels=labels)

#write.table(cpm(DG), file = paste0(prefix, "_original_cpm"), quote = F, sep = "\t")
#write.table(cpm(DG, log = T, prior_count = 1), file = paste0(prefix, "_original_logcpm"), quote = F, sep = "\t")



# Gene ID Conversion #

# Gene ID conversion using gtf used in mapping. The total number and order of genes in the gtf should be exactly the same as in the individual count files.
# Thus, can probably directly replace gene ids in DG with list of gene names in the gtf.
if (grepl(pattern = "ENSG", rownames(DG$counts)[1])) {
	gtf <- read.table(gtf_path, sep = " ")
	gene_names <- as.character(gtf[, 6])
	gene_names <- substr(gene_names, start = 1, stop = nchar(gene_names) - 1) # getting rid of the ";" at the end
	# print(gene_names[1:10])
	rownames(DG) <- gene_names
}



# Filtering #

keep <- rowSums(cpm(DG) > 1) >= length(unique(donor_list)) # filtering by expression level
print(paste0("before filtering: ", nrow(DG)), quote = F)
DG <- DG[keep, ,keep.lib.sizes = FALSE]
print(paste0("after filtering: ", nrow(DG)), quote = F)


# Deduplicating # added 7/22/2021

duplicate_genes <- row.names(DG)[duplicated(row.names(DG))]
DG <- DG[!(row.names(DG) %in% duplicate_genes), ,keep.lib.sizes = F]
print(paste0(length(duplicate_genes) + length(unique(duplicate_genes)),
			 " duplicated genes removed."))


# Gene ID convresion using org.Hs.eg.db. Results in loss of genes (57905 -> 32995).
# if(grepl("ENSG", rownames(DG$counts)[1])) {
# 	symbols <- mapIds(org.Hs.eg.db, keys = rownames(DG$counts), keytype = "ENSEMBL", column="SYMBOL")
# 	rownames(DG) <- symbols
# 	symbols_to_keep <- na.omit(symbols)
# 	DG <- DG[symbols_to_keep,]
#	print(nrow(DG))
# }



# PCA plot using DESeq2 and ggplot2 #

# samples <- data.frame(labels, donor_list, conditions)
# colnames(samples) <- c("Sample", "Experiment", "Condition")
# print(samples)
#
# dds <- DESeqDataSetFromMatrix(countData = DG$counts, colData = samples, design = ~ Experiment + Condition)
# dds <- DESeq(dds)
# vsdata <- vst(dds, blind = FALSE)
#
#
# pcafilename2 <- paste0(prefix, "_PCAplot_500.png")
# png(pcafilename2, width = 800, height = 600)
# p <- plotPCA(vsdata, intgroup = "Condition") # default ntop = 500
# p <- p + geom_text(label = colnames(vsdata), vjust = 0, nudge_y = 0.3, check_overlap = T)
# print(p)
# dev.off()


## Normalization ##

DG <- calcNormFactors(DG)
DGgroups <- c(conditions) #conditions is a comma separated list of conditions i.e. c("IFN+", "IFN+", "IFN-", "IFN-")
DG <- DGEList(counts = DG, group = factor(DGgroups))
DG <- calcNormFactors(DG) # why calcNormFactors again?
Donor <- factor(c(donor_list)) #donor_list is list of the donor of each file;
Condition <- factor(c(conditions)) #conditions is same as conditions for DGgroups;
#labels <- c(file_labels); #file_labels is what the name of each file should be when plotted;


## MDS Plot ##

mdsfilename <- paste0(prefix, "_MDSplot.pdf")
pdf(mdsfilename)

#factor_labels <- levels(factor(labels))
#print(factor_labels)

if(length(labels) > 2) {
  colors <- brewer.pal(length(labels),"Set1")
} else {
  colors <- c('deepskyblue', 'red')
}

#print(colors)

names(colors) <- labels
plotMDS(DG, col = colors, pch = 18, cex = 2)
legend('topright', title="Experiment", legend=labels,
	   col=colors, cex=.75, pch=18, xpd=T, bty="n", inset=c(0, -.162))

dev.off()


## Dispersion Estimation, BCV Plot, and Exact Test ##

data.frame(Sample=labels, Donor, Condition) # what is this line for?
design <- model.matrix(~Donor+Condition)
DG <- estimateDisp(DG, design, robust=TRUE)


bcv_file_name <- paste0(prefix, "_BCVplot.pdf")
pdf(bcv_file_name) #name of file to write bcv graph to;
plotBCV(DG) #!labels appears to be not working!;
dev.off()


#de <- exactTest(DG, pair=rev(levels(Condition))) #Only applies for one condition for two states like + and - IFN, need to use qlf or glm for others;
de <- exactTest(DG, pair = head(conditions, 2)) # the order of fold change calculation is determined by the order of first 2 variables in "variable" in edgeR_submission.sh, and the reference is the first variable.
de.genes <- rownames(topTags(de, n = n)$table)
diffex_file <- paste0(prefix, "_diffexgenes")
diffex_file_full <- paste0(prefix, "_diffexgenes_full")
write.table(de.genes, diffex_file, quote=FALSE, row.names=FALSE, col.names=TRUE) #file to write de.genes to
write.table(topTags(de, n=n)$table, diffex_file_full,
			quote=FALSE, row.names=TRUE, col.names=TRUE, sep='\t') # why do diffex_genes and diffex_genes_full have all genes instead of just n genes?
diffextable <- topTags(de, n=n)$table


## Smear Plot ##

smearplotname <- paste0(prefix, "_SmearPlot.pdf")
pdf(smearplotname)
plotSmear(DG, de.tags=de.genes) #Generates a smear plot for differentially expressed genes, upregulated and downregulated;
dev.off()

print('>>>>>>>>>>>>')


## Volcano Plot ##

detags <- topTags(de, n = n) #Choose n to be how many genes you want to show on the graph
print(colnames(detags$table))

volcanodata <- cbind(detags$table$logFC, -log10(detags$table$PValue), detags$table$FDR)
row.names(volcanodata) <- row.names(detags$table)
colnames(volcanodata) <- c("logFC", "PValue", "FDR")
write.table(volcanodata, file = paste0(prefix, "_volcanodata"),
			quote = F, row.names = T, col.names = T, sep = "\t")


volcanodata_insig <- volcanodata[volcanodata[,2] < -log10(.05),]
volcanodata_sig <- volcanodata[volcanodata[,2] >= -log10(.05),]
volcano_logFC <- volcanodata_sig[abs(volcanodata_sig[,1]) >= 2,]
volcanoplotname <- paste0(prefix, "_VolcanoPlot.png")


windowsize <- 15
#upperlim=220
print('>>>>>>>>>>>>')
print(max(volcanodata[is.finite(volcanodata[,2]),2]))
upperlim <- max(volcanodata[is.finite(volcanodata[,2]),2]) + 1
#voltitle <- paste('DEG (',length(rownames(volcano_logFC)),')')
#upreg <- paste('UpReg: ',)
#downreg <- paste('DownReg: ',)

png(volcanoplotname, height=1000, width=1000)
par(mar = c(5,5,5,5))
plot(volcanodata_insig, col='darkgray', pch=19, cex=1.5,
	 ylim=c(0,upperlim), xlim=c(-windowsize, windowsize),
	 cex.lab=2, cex.axis=1.5, font.lab=2, font=2) #Generates a volcano plot of the DE data;
points(volcanodata_sig, col='firebrick', cex=1.5, pch=19)
points(volcano_logFC, col='dodgerblue', cex=1.5, pch=19)
#legend(windowsize*0.7, upperlim*0.7, legend=c(upreg,downreg), col=c("firebrick", "dodgerblue"), pch = c(19,19), cex=1.2, title=voltitle)

#text(volcanodata_sig[volcanodata_sig[,1] > 1 & volcanodata_sig[,3] > .1, ], labels=row.names(volcanodata_sig[volcanodata_sig[,1] > 1 & volcanodata_sig[,3] > .1,]), cex=1.5, pos=1, font=2)
#Uncomment above line for COVID

#print(dim(volcanodata))
sig_FDR <- volcanodata[volcanodata[,3] < .05 & volcanodata[,2] < upperlim, ]
sig_FDR <- sig_FDR[order(abs(sig_FDR[,1]), decreasing=TRUE), ]
sig_FDR <- sig_FDR[sig_FDR[,1] < windowsize, ]
outlier <- volcanodata[volcanodata[,2] > 10, ] # <<< Adjust this number depending on the number of DE genes with higher -log10(P Value)
#print(outlier)
print(dim(outlier))
#class(outlier)
# write.table(outlier, file = paste0(prefix, "_outlier"), quote = F, row.names = T, col.names = T, sep = "\t")

num_gene <- min(15, length(rownames(sig_FDR)))
num_out <- min(20, length(rownames(outlier)))

text(sig_FDR[1:num_gene,1], sig_FDR[1:num_gene,2] ,
	 labels=row.names(sig_FDR[1:num_gene,]),
	 cex=1.2, pos=2, font=2) #Change cex to 2.5 for covid
#text(outlier[,1], outlier[,2] , labels=row.names(outlier), cex=1.2, pos=2, font=2) #Change cex to 2.5 for covid
text(outlier[1:num_out,1], outlier[1:num_out,2] ,
	 labels=row.names(outlier[1:num_out,]),
	 cex=1.2, pos=2, font=2) #Change cex to 2.5 for covid
abline(h=-log10(.05), col="blue", lty=2, lwd=3)

dev.off()


## Heatmap ##

y <- cpm(DG, log=TRUE, prior.count = 1)
cpmname <- paste0(prefix, "_cpm")
write.table(y, cpmname, quote=FALSE, row.names=TRUE, col.names=TRUE, sep = "\t")
volcdata_sub <- volcanodata[rownames(volcanodata) %in% rownames(y),]
max_rows <- rownames(volcdata_sub[order(abs(volcdata_sub[volcdata_sub[,3] < .05, ][,1]), decreasing=TRUE), ]) # filter genes with FDR <0.05, and order according to logFC.
head(max_rows)
#selY <- y[rownames(tbl)[tbl$FDR<.01 & abs(tbl$logFC) > 1.5],]; #Changed significant from .1 to .01
selY <- y[max_rows, ]
selY <- selY[1:min(30, length(rownames(selY))),] # selects top 30 genes with the most logFC.

htmaptablename <- paste0(prefix, "_heatmap_table")
write.table(selY, htmaptablename, quote=FALSE, sep = "\t")

heatmapname <- paste0(prefix, "_HeatMap.png")
png(heatmapname, height=650, width=500) #Generates a heatmap showing how expression changes of top DE genes;
#heatmap.2(selY, col=bluered(100), , trace="none", density.info="none", symkey=F, scale="none", cexRow=1.5,cexCol=2,margins=c(14,16), Colv=FALSE);
Heatmap(selY,
		column_order=colnames(selY),
		heatmap_height=unit(20,'cm'),
		heatmap_width=unit(12,'cm'),
		heatmap_legend_param = list(title = "log2(CPM)"))
dev.off()


## Collapsed Heatmap ##

cols <- unique(conditions) # order of conditions in collapsed heatmap depends on order in conditions, or mock first.
collapsed <- matrix(nrow=nrow(selY), ncol=length(cols))
if("mock" %in% cols) {
	cols <- c("mock", cols[cols != "mock"])
}
colnames(collapsed) <- cols
cols <- colnames(collapsed)
rownames(collapsed) <- rownames(selY)

for(i in seq(1, length(cols))) {
	collapsed[,i] <- rowMeans(selY[,grepl(cols[i], colnames(selY))])
}

collapsename <- paste(prefix, "_collapsed_heatmap.png", sep="")
png(collapsename, height=650, width=350) #Generates a heatmap showing how expression changes of top DE genes;
#heatmap.2(collapsed, col=bluered(100), trace="none", lhei=c(1,10), lwid=c(1,2), srtCol=0, adjCol=c(0.5,0), offsetCol=2, density.info="none", symkey=F, scale="none", cexRow=1.5,cexCol=2,margins=c(10,10), Colv=FALSE);
Heatmap(collapsed,
		column_order=colnames(collapsed),
		heatmap_height=unit(20,'cm'),
		heatmap_width=unit(8,'cm'),
		heatmap_legend_param = list(title = "log2(CPM)"))
dev.off()


# mock_vs_exp_df <- data.frame(mock=.5*(selY[,3] + selY[,1]), IFIT2=.5*(selY[,4] - selY[,2]))

# mock_v_exp <- paste(prefix, "_mock_IFIT_HeatMap.png", sep="")
# png(mock_v_exp, height=1000, width=1000); #Generates a heatmap showing how expression changes of top DE genes;
# heatmap.2(as.matrix(mock_vs_exp_df), col=redgreen(75), , trace="none", density.info="none", symkey=F, scale="none", cexRow=2,cexCol=2,margins=c(14,16), Colv=FALSE);
# dev.off();


## LogFC Heatmap ##

#logFC_df <- data.frame(exp32_logFC=selY[,2] - selY[,1], exp37_logFC=selY[,4] - selY[,3])

donors <- unique(donor_list)
logFC_df <- matrix(nrow=nrow(selY), ncol=length(donors), dimnames=list(rownames(selY), donors))
for (i in seq(length(donors))) {
	logFC_df[,i] <- selY[,2*i] - selY[,2*i-1] # logFC calculation depends the order of variables (columns) in selY, which depends on filelist order.
}

logFC_name <- paste(prefix, "_logFC_HeatMap.png", sep="")
png(logFC_name, height=650, width=350) #Generates a heatmap showing how expression changes of top DE genes;
#heatmap.2(as.matrix(logFC_df), col=bluered(100), , trace="none", density.info="none", symkey=F, scale="none", cexRow=2,cexCol=2,margins=c(14,16), Colv=FALSE);
#heatmap.2(diffextable[1:80], col=bluered(100), , trace="none", density.info="none", symkey=F, scale="none", cexRow=2,cexCol=2,margins=c(14,16), Colv=FALSE);
#Heatmap(as.matrix(logFC_df), column_order=colnames(logFC_df), heatmap_height=unit(20,'cm'), heatmap_width=unit(8,'cm'))
Heatmap(logFC_df,
		column_order=colnames(logFC_df),
		heatmap_height=unit(20,'cm'),
		heatmap_width=unit(8,'cm'),
		heatmap_legend_param = list(title = "log2(CPM)"))
dev.off()


# heatmap with specific geneset (ISGs) ##

# isgs <- read.table('/home/skutluay/Test_Files/nakys_things/Scripts/isgs.txt',header=F)
# gset <- intersect(as.character(isgs[,1]),rownames(y))
# print('///////////////////////')
# print(dim(y[gset,]))
# #print(class(y[gset,]))
# geneset_plot_name <- paste(prefix, "_target_geneset.png", sep="")
# png(geneset_plot_name, height=650, width=500)
# #heatmap.2(y[gset,], col=bluered(100), lhei=c(1,10), lwid=c(1,2), trace="none", density.info="none", symkey=F, scale="none", cexRow=1.5,cexCol=2,margins=c(14,16), Colv=FALSE)
# Heatmap(y[gset,], column_order=colnames(y), heatmap_height=unit(20,'cm'), heatmap_width=unit(12,'cm'))
# dev.off()
