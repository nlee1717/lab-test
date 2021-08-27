# Functions for ribosome profiling analysis

make_output_directories <- function(outputDir) {
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "objs"))) {
    dir.create(file.path(outputDir, "objs"), recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "reports", "figs"))) {
    dir.create(file.path(outputDir, "reports", "figs"), recursive = T)
  }
  if (!dir.exists(file.path(outputDir, "reports", "de_genes"))) {
    dir.create(file.path(outputDir, "reports", "de_genes"), recursive = T)
  }
  
}


make_dge <- function(counts_path, countFilter, region="cds") {
  # function to make dge object for Ribo-seq or RNA-seq and output the counts table and CPM table
  # Args: 
  #   counts_path: path for the output count tables from htseq-count
  #   countFilter: common pattern of the htseq file names in the data type
  #   region: cds, utr3, or utr5
  # Returns:
  #   null. But will save the dge object 
  
  files <- dir(path = counts_path, pattern = "*.count_reduced$")
  
  ## sort the files in asc order
  files <- sort(files)
  counts <- readDGE(files, path = counts_path, header = F)$counts
  print(dim(counts))
  # filter out metadata fields in counts (not currently useful)
  noint <- rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual",
                                   "__not_aligned","__alignment_not_unique")
  counts <- counts[!noint,]
  print(noquote("counts:"))
  print(dim(counts))
  #print(colnames(counts))
  
  # get gene information from ensembl
  #ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", version = "Ensembl Genes 75", dataset = "hsapiens_gene_ensembl", host = "http://feb2014.archive.ensembl.org") # human
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", version = "Ensembl Genes 90", dataset = "csabaeus_gene_ensembl", host = "http://aug2017.archive.ensembl.org") # Vero
  
  stopifnot(length(unique(rownames(counts))) == nrow(counts))
  
  # gtf.ens <- getBM(attributes=c('ensembl_gene_id','external_gene_id',"gene_biotype"), filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl) # human
  gtf.ens <- getBM(attributes=c('ensembl_gene_id','external_gene_name',"gene_biotype"), filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl) # Vero
  gtf.ens <- rename(gtf.ens, external_gene_id = external_gene_name) # Vero
  
  saveRDS(gtf.ens, file.path(file.path(OUTPUT, "objs"), paste0(region, "_genes.rds")))
  
  genes <- readRDS(file.path(file.path(OUTPUT, "objs"), paste0(region, "_genes.rds")))
  print(noquote("genes:"))
  print(dim(genes))
  
  genes <- genes[match(rownames(counts), genes$ensembl_gene_id),]

  # added lines ##############################
  genes <- na.omit(genes)
  counts <- counts[genes$ensembl_gene_id, ] # human
  ############################################
  
  stopifnot(genes$ensembl_gene_id == rownames(counts))

  ## get subset of columns from count table
  counts_sub <- counts[, grep(countFilter, colnames(counts))]

  ## get samples information
  samples <- read_csv("samples.csv")
  print(samples)
  # make sure order of the colnames in counts is the same as samples
  counts_sub <- counts_sub[, paste0(samples$sample, "")]
  #stopifnot(str_remove(colnames(counts_sub), "") == samples$sample)
  stopifnot(colnames(counts_sub) == samples$sample)
  write.csv(counts_sub, file.path(file.path(OUTPUT, "reports"), paste0("gene_counts_", region, ".csv")))
  colnames(counts_sub) <- paste(samples$experiment, samples$independent_var, samples$experiment_type, sep = "_")
  dge <- DGEList(counts = counts_sub, 
                 group = paste("time", samples$independent_var, samples$experiment_type, sep = "_"),
                 genes = genes)
  dge$samples$experiment <- samples$experiment
  dge$samples$independent_var <- samples$independent_var
  dge$samples$experiment_type <- samples$experiment_type
  
  # store only protein coding genes for later analysis
  dge <- dge[dge$genes$gene_biotype == "protein_coding",]
  dge <- na.omit(dge)
  print(noquote("protein coding: "))
  print(dim(dge$genes))
  
  dge <- calcNormFactors(dge)
  saveRDS(dge, file.path(file.path(OUTPUT, "objs"), paste0("dge_", region, "_protein_coding.rds")))
  
  cpm <- cpm(dge, log = T, prior.count = 1)
  cpm <- as.data.frame(cpm)
  cpm$ensembl_gene_id <- rownames(cpm)
  
  write.csv(as.data.frame(cpm), file.path(file.path(OUTPUT, "reports"), paste0("cpm_", region, ".csv")), row.names = F)
  write.csv(colSums(dge$counts), file.path(file.path(OUTPUT, "reports"), paste0("colsums", region, ".csv")))
}

generate_reprod_plot <- function(cpmTable, dataType) {
  # function to generate reproducibility plot
  # Args:
  # cpmTable: cpm table created by make_dge
  # dataType: experiment name either rnaseq or riboseq
  
  cor_ <- function(df, col1, col2) {
    # if this line returns error, see if df colnames is infected or timepoint
    stopifnot(c(col1, col2) %in% names(df))
    return(paste0("r = ", format(cor(df[[col1]], df[[col2]]), digits = 3)))
  }
  colors<-list("rnaseq" = "deepskyblue", "riboseq" = "darkorchid2")
  samples <- paste(SAMPLES, dataType, sep = "_")
  comb <- matrix(samples, nrow=2)
  plots <- list()
  for(i in 1:ncol(comb)) { 
    df.cor <- cor_(cpmTable, comb[1, i], comb[2, i])
    plots[[i]] <- ggplot(cpmTable, aes_string(paste0('`', comb[1, i], '`'), 
                                              paste0('`', comb[2, i], '`'))) +
      geom_point(alpha=.1, colour = colors[[dataType]]) +
      annotate("text", x=9, y = 0, label = df.cor, size = 4) +
      labs(x=paste(dataType, comb[1, i], sep = " "), y=paste(dataType, comb[2, i], sep = " ")) +
      theme_bw(base_size = 8) 
  }
  if (save_figs) pdf(file.path(OUTPUT, "reports", "figs", paste0(dataType, "_reproducibility.pdf")))
  do.call(grid.arrange, c(plots, ncol=4))
  if (save_figs) dev.off()
}

generate_mds <- function(dge, dataType) {
  # function to generate MDS plot
  # Args:
  #   dge: DGE object create by make_dge
  #   dataType: experiment name either rnaseq or riboseq
  
  colors<-c('deepskyblue', 'darkorchid2', 'seagreen2', 'midnightblue')
  dge_sub <- dge[, grep(dataType, colnames(dge))]
  if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(dataType, "_MDS_dim_1_2.pdf")))
  plotMDS(dge_sub, col=colors, dim.plot = c(1,2))
  #legend('topright', legend=dge_sub$samples$group, col=colors, cex=1, pch=18)
  if (save_figs) dev.off()
  
  if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(dataType, "_MDS_dim_1_3.pdf")))
  plotMDS(dge_sub, col=colors, dim.plot = c(1,3))
  #legend('bottomright', legend=dge_sub$samples$group, col=colors, cex=1, pch=18)
  if (save_figs) dev.off()
  
  if (save_figs) pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(dataType, "_MDS_dim_2_3.pdf")))
  plotMDS(dge_sub, col=colors, dim.plot = c(2, 3))
  #legend('bottomright', legend=dge_sub$samples$group, col=colors, cex=1, pch=18)
  if (save_figs) dev.off()
}

plot_CV_vs_CPM <- function(cpmTable, repList, dataType) {
  # function to plot the coefficient variation vs the mean log2 CPM of replications 
  # Args:
  #   cpmTable: CPM table generated by make_dge
  #   repList: specify the list of replications
  #   dataType: either rnaseq or riboseq
  for (reps in repList) {
    reps <- paste(reps, dataType, sep = "_")
    tb <- cpmTable %>% 
      mutate(mean=abs(rowMeans(.[reps])), sd = rowSds(.[reps]), cv = sd/mean,
             reps = paste(reps, collapse = "_vs_"))
    reps_name <- paste(reps, collapse = "_and_")
    
    pt <- ggplot(tb, aes(x=mean, y=cv)) +
      geom_point(size=0.5) +
      ylim(c(0, 25)) +
      xlim(c(0, 25)) +
      labs(x = "Mean log2 CPM", y = "Coefficient variantion", title = reps_name)
    
    if (save_figs) {
      pdf(file.path(file.path(OUTPUT, "reports", "figs"), paste0(reps_name, "_CV_vs_CPM.pdf")))
      print(pt)
      dev.off()
    }
    
    print(pt)
  }
  
}

filter_dge <- function(dge, region, ribothreshold = 1, rseqthreshold = 1) {
  # function to filter out genes with very low counts across all libraries
  # Genes must have at least NUMREP number of ribo-seq samples with count more than
  # RIBOTHRESHOLD, and at least NUMREP number of rna-seq samples with count more than
  # RSEQTHRESHOLD, to pass the filter
  # Args
  #   dge: DGE obejct created by make_dge
  #   region: cds, utr3 or utr5
  #   ribothreshold: cut-off cpm value for riboseq samples for filtering
  #   rseqthreshold: cut-off cpm value for rnaseq samples for filtering
  #   #threshold: cut-off value of log2 CPM for filtering
  #   #numRep: minimum number of replicates
  print("dimension before filtering:")
  print(dim(dge))
  #print(head(cpm(dge, log=T, prior.count = 1), n=10))
  #write.csv(cpm(dge), file = file.path(OUTPUT, "reports", paste0("original_cpm_", region, ".csv")))
  
  numrep <- ncol(dge) / 4
  #keep <- rowSums(cpm(dge) > threshold) >= numRep 
  cpm_temp <- cpm(dge)
  keep <- (rowSums(cpm_temp[, 1:(ncol(cpm_temp) / 2)] > ribothreshold) >= numrep) &
    (rowSums(cpm_temp[, (ncol(cpm_temp) / 2 + 1):(ncol(cpm_temp))] > rseqthreshold) >= numrep)
  
  print(table(keep))
  dge <- dge[keep, ,keep.lib.sizes=FALSE]
  print("dimension after filtering:")
  print(dim(dge))
  #print(head(dge))
  
  # normalization
  dge <- calcNormFactors(dge)
  saveRDS(dge, file.path(OUTPUT, "objs", paste0("dge_", region, "_filt_norm.rds")))
  cpm <- cpm(dge, log = T, prior.count = 1) # default prior.count is 2
  cpm <- as.data.frame(cpm)
  cpm$ensembl_gene_id <- rownames(cpm)
  write.csv(cpm, file.path(file.path(OUTPUT, "reports"), paste0("cpm_", region, "_filt.csv")), row.names = F)
}

get_table_for_samples <- function(samplesList, cpmTable) {
  # function to calculate TE of each sample 
  # Args:
  #   samplesList: names of all samples without dataType 
  #   cpmTable: CPM table from make_dge
  
  table <- data.frame(ensembl_gene_id = cpmTable$ensembl_gene_id)
  
  for(sample in samplesList) {
    df <- data.frame(ensembl_gene_id = cpmTable$ensembl_gene_id, rnaseq = cpmTable[[paste0(sample, "_rnaseq")]], rf = cpmTable[[paste0(sample, "_riboseq")]])
    table <- table %>%
      left_join(df, by = "ensembl_gene_id") %>%
      #mutate(TE = rf / rnaseq)
      mutate(TE = rf - rnaseq)
    colnames(table)[colnames(table) == "rnaseq"] <- paste0(sample, "_rnaseq")
    colnames(table)[colnames(table) == "rf"] <- paste0(sample, "_rf")
    colnames(table)[colnames(table) == "TE"] <- paste0(sample, "_te")
  }
  return (table)
}

plot_individual_density <- function(combinedTable, region) {
  # function to make density plot of rnaseq, rf density plot for cds, utr3, utr5 individually
  for (sample in SAMPLES) {
    dataType <- c("rnaseq", "rf", "te")
    for (dt in dataType) {
      dens_plot <- ggplot(combinedTable, aes_string(x=paste0(sample, "_", dt))) +
        geom_density() +
        xlim(c(-20, 20)) +
        labs(title = paste0("Log2 CPM density for sample ", dt, " in region ", region))
      print(dens_plot)
    }
  }
}

de_analysis <- function(dge, region = "") {
  # function to do different expression analysis with DGE object and output de genes to reports/de_genes, using edgeR.
  # Args:
  #   dge: DGE object from make_dge
  ## Make contrast
  
  
  group <- factor(dge$samples$group)
  print(group)
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  contrasts <- makeContrasts(
    time_infected_vs_mock_rnaseq = time_infected_rnaseq - time_mock_rnaseq,
    time_infected_vs_mock_riboseq = time_infected_riboseq - time_mock_riboseq,
    time_infected_vs_mock_te = (time_infected_riboseq - time_mock_riboseq) - (time_infected_rnaseq - time_mock_rnaseq),
    levels = design
  )
  
  
  dge <- estimateDisp(dge, design)
  ##  A common dispersion (i.e. red line on the BCV plot) between 0.2 and 0.4 is usually considered reasonable and hence could detect more DE genes.If the common dispersion is above the 0.4 threshold, this will influence the number of DE genes found in the study. 
  # note: red line on BCV plot is biological coefficient of variation (square root of common dispersion), not common dispersion.
  message("common dispersion = ", dge$common.dispersion)
  plotBCV(dge)
  
  fit <- glmFit(dge, design)
  
  getDE <- function(contrast_name) {
    topTags(glmLRT(fit, contrast=contrasts[, contrast_name]), n=Inf, sort.by = "PValue")[[1]]
  }
  de <- lapply(attr(contrasts, 'dimnames')$Contrasts, getDE)
  
  names(de) <- attr(contrasts, 'dimnames')$Contrasts
  saveRDS(de, file.path(OUTPUT, "objs", paste0("de", region, ".rds")))
  
  de_results_path <- file.path(OUTPUT, "reports", "de_genes")
  for (comp in attr(contrasts, 'dimnames')$Contrasts) {
    write_csv(filter(de[[comp]], FDR < .1), file.path(de_results_path, paste0(comp, region, ".csv")))
    write_csv(de[[comp]], file.path(de_results_path, paste0(comp, region, "_full.csv")))
  }
}

plot_heatmap <- function(table, save_fig = F, figname = "heatmap") {
  # function to plot heatmap of log fold change of de genes
  # Args:
  #   table: put log FC of de genes from each comparison into one table. The genes are the union of all de genes

  
  mat <- as.matrix(table[, -grep("ensembl_gene_id", colnames(table))])
  # set NA to zero
  mat[is.na(mat)] <- 0
  breaks <- c(-6, -4, -2, -1, 0, 1, 2, 4, 6)
  mypalette <- rev(brewer.pal(8, "RdBu"))
  if (save_fig) pdf(file.path(OUTPUT, "reports", "figs", paste0(figname, ".pdf")))
  heatmap.2(mat,
            dendrogram = "none",
            trace = "none",
            density.info = "none",
            labRow = "", 
            key.xlab = "Log2 FC",
            key.title = NA,
            margins = c(10, 5),
            col = mypalette,
            cexCol = 1,
            breaks = breaks,
            symkey = F,
            Colv = F
  )
  if (save_fig) dev.off()
}


