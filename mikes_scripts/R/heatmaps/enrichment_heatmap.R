# Gene set over-representation analysis for gene clusters from time course cpm heatmap

library(dplyr)
library(msigdbr)
library(clusterProfiler)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)

basepath <- "C:/kutluaylab/data/Kyle_Exp92/Ribo-seq"
path <- file.path(basepath, "enrichment_heatmap")
cpm_path <- file.path(basepath, "heatmap", "sig_zscore")
prefix <- "ribo_all_heatmap"


### DATA PREPARATION ###

# filtered cpm table from time course cpm heatmap for gene universe
cpm <- read.delim(file.path(cpm_path, paste0(prefix, "_cpm")),
                  row.names = NULL)

prefix <- "ribo_all_heatmap_v2"
outputname <- paste0(prefix, "_ER")

# cpm table of DE genes for cluster information
sig_cpm <- read.delim(file.path(cpm_path, paste0(prefix, "_sig_cpm_cluster")))



gene_cluster <- na.omit(data.frame(row.names(sig_cpm), sig_cpm$cluster))
colnames(gene_cluster) <- c("gene_name", "cluster")
# total number of clusters
num_cluster <- max(gene_cluster$cluster)
selected_clusters <- seq(num_cluster)


cpm <- rename(cpm, gene_name = row.names)
duplicate_genes <- cpm$gene_name[duplicated(cpm$gene_name)] # getting rid of duplicate genes in cpm table, including the original one
cpm <- cpm[!(cpm$gene_name %in% duplicate_genes), ]
print(paste0(length(duplicate_genes) + length(unique(duplicate_genes)),
             " duplicated genes removed."))



geneNames.universe <- cpm$gene_name


# getting rid of genes in the duplicate list that may be present in clustered genes
gene_cluster <- gene_cluster[!gene_cluster$gene_name %in% duplicate_genes, ]
stopifnot(gene_cluster$gene_name %in% geneNames.universe)

### just for short Vero riboseq ###
selected_clusters <- c(1, 2, 4, 5, 6)
gene_cluster_orig <- gene_cluster
gene_cluster <- gene_cluster[gene_cluster$cluster %in% selected_clusters, ]
num_cluster <- length(selected_clusters)
###################################


## Retrieving info from Molecular Signatures database ##
msig_bp <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H" | (gs_cat == "C5" & gs_subcat == "GO:BP")) %>%
  select(gs_name, human_gene_symbol)

save(msig, file = file.path(path, "enrichment_results", "ER_heatmap_msig.RData"))
load(file = file.path(path, "enrichment_results", "ER_heatmap_msig.RData"))

msig_mf <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "C5" & gs_subcat == "GO:MF") %>%
  select(gs_name, human_gene_symbol)

msig_tft <- msigdbr() %>%
  filter(gs_cat == "C3" & gs_subcat == "TFT:GTRD") %>%
  select(gs_name, human_gene_symbol)


msig <- msigdbr() %>%
  filter(gs_cat == "H" |
           (gs_cat == "C5" & gs_subcat == "GO:BP") |
           (gs_cat == "C5" & gs_subcat == "GO:MF") |
           (gs_cat == "C3" & gs_subcat == "TFT:GTRD")) %>%
  select(gs_name, human_gene_symbol)

save(msig, file = file.path(path, "enrichment_results", "ER_heatmap_msig.RData"))
load(file = file.path(path, "enrichment_results", "ER_heatmap_msig.RData"))




### OVER-REPRESENTATION ANALYSIS ###

enrich <- function(cluster_table, universe, database, selected_clusters) {
  enrichment.results <- NULL

  ## Over-representation analysis for each cluster, filtering q-value <= 0.1 ##
  for (cluster.id in selected_clusters) {

    geneNames.cluster <- cluster_table[cluster_table$cluster == cluster.id, ]$gene_name

    em <- enricher(gene = geneNames.cluster,
                   pAdjustMethod = "fdr",
                   universe = universe,
                   TERM2GENE = database,
                   minGSSize = 10,
                   qvalueCutoff = 0.1)


    if (!is.null(em)) {
      result <- em@result
      result$cluster <- cluster.id
      result <- filter(result, qvalue <= 0.1)
      if (nrow(result) > 0) {
        # result$genes <- apply(result, 1, function(x) {str_c(sort(getSYMBOL(unlist(str_split(x[["geneID"]],"/")), data='org.Hs.eg')),collapse=",")} )
        enrichment.results <- rbind(enrichment.results, result)
      }
    }
  }

  return(enrichment.results)
}


enrichment.results <- enrich(cluster_table = gene_cluster,
                             universe = geneNames.universe,
                             database = msig,
                             selected_clusters = selected_clusters)


# enrichment.results_mf <- enrich(gene_cluster = gene_cluster,
#                                 geneNames.universe = geneNames.universe,
#                                 database = msig_mf,
#                                 selected_clusters = selected_clusters)
#
# enrichment.results_tft <- enrich(gene_cluster = gene_cluster,
#                                  geneNames.universe = geneNames.universe,
#                                  database = msig_tft,
#                                  selected_clusters = selected_clusters)




## Function to run enricher for individual clusters ##
enrichCluster <- function(gene_cluster, cluster, universe, msig, path) {
  genes <- gene_cluster[gene_cluster$cluster == cluster, ]$gene_name
  enrichment <- enricher(gene = genes,
                         pAdjustMethod = "fdr",
                         universe = universe,
                         TERM2GENE = msig,
                         minGSSize = 10,
                         qvalueCutoff = 0.1)

 enrichment_result <- enrichment@result
 write.csv(enrichment_result, file.path(path, "enrichment_results", paste0("cluster", cluster, "_results.csv")))
}

enrichment.results$cluster <- as.integer(enrichment.results$cluster)
enriched_cluster <- unique(enrichment.results$cluster)

enrichment.results <- rename(enrichment.results, genes = geneID)

enrichment.results$category <- ""
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^HALLMARK")] <- "Hallmark"
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^GOBP")] <- "Biological Process"
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "^GOMF")] <- "Molecular Function"
enrichment.results$category[str_detect(enrichment.results$Description, pattern = "TARGET")] <- "Transcription Factor"



## Getting number of genes present in each cluster ##
num.in.cluster <- data.frame(table(gene_cluster$cluster))
colnames(num.in.cluster) <- c("cluster", "n")
write.csv(num.in.cluster,
          file.path(path, "enrichment_results", paste0(outputname, "_num_in_cluster.csv")),
          row.names = F)


## Getting the ratio of genes in a cluster that belong to a term against the total number of genes in the cluster ##
enrichment.results$ratio <- enrichment.results$Count / num.in.cluster[match(enrichment.results$cluster, num.in.cluster$cluster), ]$n
write.csv(enrichment.results,
          file.path(path, "enrichment_results", paste0(outputname, "_results.csv")),
          row.names = F)


## "Flattening" the results by showing results for each cluster in parallel, merging terms shared by different clusters ##
enrichment.results_wide <- pivot_wider(enrichment.results,
                                       id_cols = c(Description, category),
                                       names_from = cluster,
                                       values_from = c(pvalue, p.adjust, qvalue, Count,
                                                       GeneRatio, BgRatio, ratio, genes))
enrichment.results_wide <- arrange(enrichment.results_wide, category)
write.csv(enrichment.results_wide,
          file.path(path, "enrichment_results", paste0(outputname, "_results_wide.csv")))
save(enrichment.results_wide,
     file = file.path(path, "enrichment_results", paste0(outputname, "_results_wide.RData")))

# Kyle's rseq v1 stopped here #



## Manual selection of interesting terms ##

outputname <- "ribo_all_heatmap_v2_ER"
load(file.path(path, "enrichment_results", paste0(outputname, "_results_wide.RData")))


selected_term_index <- read.table(file = file.path(path, "selected_term", "selected_term_index.txt"),
                                  header = F)[, 1] # selected manually

enrichment.results.pruned <- enrichment.results_wide[selected_term_index, ]

### Vero ribo short only ###
selected_term_upper <- read.table(file = file.path(path, "selected_term", "selected_term_upper.txt"),
                                  header = F)[, 1]
enrichment.results.pruned <- enrichment.results_wide[enrichment.results_wide$Description %in% selected_term_upper, ]
excluded_terms <- setdiff(selected_term_upper, enrichment.results.pruned$Description)
######################

enrichment.results.pruned <- arrange(enrichment.results.pruned, category)

# clusters that have selected terms
selected_clusters <- c(1, 4)
selected_clusters <- enriched_cluster
## for Kyle Exp92 ribo only, to only include clusters with selected terms, need improvement ##
enrichment.results.pruned <- enrichment.results.pruned[, c(1, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17,
                                                           18, 20, 21, 23, 24, 26)]
###################################


## Changing the style of the Description column ##
selected_terms <- tolower(enrichment.results.pruned$Description) %>%
  gsub(pattern = "_", replacement = " ") %>%
  substr(start = nchar(word(., 1)) + 2, stop = nchar(.))

write.table(data.frame(selected_terms),
            file = file.path(path, "selected_term", paste0(outputname, "_selected_terms.txt")),
            col.names = F, row.names = F, quote = F)


## After manual editing of descriptions ##
selected_terms_edited <- read.table(file = file.path(path, "selected_term", paste0(outputname, "_selected_terms_edited.txt")),
                                    header = F,
                                    sep = "\t")[, 1]

enrichment.results.pruned$Description <- selected_terms_edited


## Only for HBEC ribo ##
enrichment.results.pruned[22, "Description"] <- "Coagulation (HM)"
enrichment.results.pruned[22, "category"] <- "Biological Process"
enrichment.results.pruned[23, "Description"] <- "KRAS Signaling Up (HM)"
enrichment.results.pruned[23, "category"] <- "Biological Process"

custom.order <- c(1, 2, 3, 4, 5, 6, 7, 8, 22, 23, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
enrichment.results.pruned <- enrichment.results.pruned[custom.order, ]
########################

## Only for Kyle Exp92 rseq ##
enrichment.results.pruned[19, "Description"] <- "DNA Binding Transcription Activator Activity (MF)"
enrichment.results.pruned[19, "category"] <- "Biological Process"
enrichment.results.pruned[20, "Description"] <- "DNA Binding Transcription Repressor Activity (MF)"
enrichment.results.pruned[20, "category"] <- "Biological Process"

custom.order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 20, 16, 17, 18)
enrichment.results.pruned <- enrichment.results.pruned[custom.order, ]
#############################

## Only for Kyle Exp92 ribo ##
enrichment.results.pruned[12, "Description"] <- paste0(enrichment.results.pruned[12, "Description"], " (HM)")
enrichment.results.pruned[13, "Description"] <- paste0(enrichment.results.pruned[13, "Description"], " (MF)")
enrichment.results.pruned[14, "Description"] <- paste0(enrichment.results.pruned[14, "Description"], " (MF)")
enrichment.results.pruned[15, "Description"] <- paste0(enrichment.results.pruned[15, "Description"], " (MF)")
# did not change those rows to "Biological Process" here
##############################


save(enrichment.results.pruned,
     file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.RData")))
load(file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.RData")))

write.csv(enrichment.results.pruned,
          file = file.path(path, "enrichment_results", paste0(outputname, "_pruned.csv")),
          row.names = F)



## Selecting q-values for plotting ##
em.q.values <- select(enrichment.results.pruned, grep(pattern = "qvalue",
                                                      colnames(enrichment.results.pruned),
                                                      value = T))
# enrichment.results.pruned$Description[20] <- "Inflammatory Response " # Vero rseq

enriched_cluster <- 1:5 # HBEC rseq only
em.q.values <- em.q.values[, 1:5] # HBEC rseq only, to get rid of cluster 6 which has no selected terms

em.q.values <- as.data.frame(em.q.values)
row.names(em.q.values) <- as.character(enrichment.results.pruned$Description)
# colnames(em.q.values) <- seq(num_cluster)
colnames(em.q.values) <- selected_clusters
  # colnames(em.q.values) <- c(2, 5) # HBEC ribo
# colnames(em.q.values) <- 1:3 # Vero rseq
# colnames(em.q.values) <- c(1, 2, 3, 5, 6) # old Vero ribo
# colnames(em.q.values) <- c(1, 2, 3, 5, 6, 7, 9, 10) # Vero ribo

em.q.values <- -log10(as.matrix(em.q.values))

write.csv(em.q.values,
          file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.csv")),
          row.names = F)

save(em.q.values,
     file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.RData")))
# load(file = file.path(path, "enrichment_results", paste0(outputname, "_qvalues.RData")))




## Selecting ratio values for plotting ##
ratio.values <- select(enrichment.results.pruned, grep(pattern = "ratio",
                                                       colnames(enrichment.results.pruned),
                                                       value = T))

ratio.values <- ratio.values[, 1:5] # HBEC rseq only
ratio.values <- as.data.frame(ratio.values)
row.names(ratio.values) <- as.character(enrichment.results.pruned$Description)
# colnames(ratio.values) <- seq(num_cluster)
# colnames(ratio.values) <- c(2, 5) # HBEC ribo
# colnames(ratio.values) <- 1:3 # Vero rseq
# colnames(ratio.values) <- c(1, 2, 3, 5, 6) # old Vero ribo
# colnames(ratio.values) <- c(1, 2, 3, 5, 6, 7, 9, 10) # Vero ribo
colnames(ratio.values) <- selected_clusters

write.csv(ratio.values,
          file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.csv")),
          row.names = F)

save(ratio.values,
     file = file.path(path, "enrichment_results", paste0(outputname, "_ratio.RData")))





### ENRICHMENT HEATMAP ###

## col_fun function for heatmap to plot colors of dots according to q-values ##
max.q <- max(em.q.values, na.rm = T)
min.q <- min(em.q.values, na.rm = T)

col_fun <- colorRamp2(c(min.q, max.q), c("black", "red")) # should depend on distribution of em.q.values data


## Column annotation blocks with cluster information ##
colPalette <- c("#AEC7E87F", "#98DF8A7F", "#FF98967F",
                "#C49C947F", "#C5B0D57F", "#FFBB787F",
                "turquoise", "pink", "maroon", "wheat3")


col.annot <- HeatmapAnnotation(cluster = anno_simple(colnames(em.q.values),
                                                     height = unit(3, "mm"),
                                                     # col = structure(colPalette[1:num_cluster],
                                                     #                 names = colnames(em.q.values))
                                                     # col = structure(colPalette[c(2, 5)],
                                                     #                 names = c(2, 5)) # HBEC ribo
                                                     # col = structure(colPalette[1:3],
                                                     #                 names = 1:3) # Vero rseq
                                                     # col = structure(colPalette[c(1, 2, 3, 5, 6)],
                                                     #                 names = c(1, 2, 3, 5, 6)) # old Vero ribo
                                                     # col = structure(colPalette[c(1, 2, 3, 5, 6, 7, 9, 10)],
                                                     #                 names = c(1, 2, 3, 5, 6, 7, 9, 10)) # Vero ribo
                                                     col = structure(colPalette[selected_clusters],
                                                                     names = selected_clusters)
                                                     ),
                               show_annotation_name = F,
                               show_legend = F)



## Heatmap ##
h <- Heatmap(em.q.values,
             # labels
             column_title = "Enriched Gene Sets",
             column_title_gp = gpar(fontsize = 9),
             column_names_gp = gpar(fontsize = 10),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize = 12),
             row_names_gp = gpar(fontsize = 9),

             # annotation
             bottom_annotation = col.annot,

             # legends
             show_heatmap_legend = F,
             # name = "-log10(q-value)",
             # col = col_fun,
             # heatmap_legend_param = list(color_bar = "continuous",
             #                             title_gp = gpar(fontsize = 12),
             #                             labels_gp = gpar(fontsize = 12),
             #                             direction = "vertical"
             #                             #legend_width = unit(4, units = "cm")
             #                             ),

             show_row_names = T,
             show_column_names = T,

             # clustering
             cluster_rows = F,
             cluster_columns = F,

             # split = enrichment.results.pruned$category,

             # circles
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height,
                         gp = gpar(fill = "white", col = "#EEEEEE")) # white rectangles in the background
               grid.circle(x = x, y = y, r = sqrt(ratio.values[i, j]) * 4, # <<< change size of circles according to ratio distribution
                           default.units = "mm",
                           # coloring the circles according to the q-values
                           gp = gpar(fill = col_fun(em.q.values[i, j]), col = NA))
             }

             # heatmap_height = unit(10, "cm"),
             # heatmap_width = unit(10, "cm")
)



## Q-value legend ##
q_legend <- Legend(col_fun = col_fun,
                   title = "-log10(q-value)",
                   title_gp = gpar(fontsize = 9),
                   labels_gp = gpar(fontsize = 8),
                   grid_height = unit(4, "mm"),
                   direction = "horizontal")



## Ratio value legend ##
ratio_breaks <- c(0.5, 0.4, 0.3, 0.2, 0.1) # <<< change manually
# ratio_breaks <- c(0.25, 0.2, 0.15, 0.1, 0.05) # Vero rseq
# ratio_breaks <- c(0.8, 0.6, 0.4, 0.2) # Vero ribo

ratio_legend <- Legend(labels = ratio_breaks,
                       grid_height = unit(6, "mm"),
                       grid_width = unit(6, "mm"),
                       title = "Percent\nof Cluster",
                       title_gp = gpar(fontsize = 9),
                       labels_gp = gpar(fontsize = 8),
                       type = "points",
                       pch = 1, # circle
                       size = unit(sqrt(ratio_breaks) * 10, "mm"), # <<< not the same as in cell_fun, change for consistency
                       # direction = "horizontal", # not working for some reason
                       # legend_gp = gpar(col = "black"),
                       background = "white"
)



# outputname <- "HBEC_ribo_ER_heatmap_3"

# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1350, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 20), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.86, "npc"))
# dev.off() # Vero rseq

# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1600, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 15), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.89, "npc"))
# dev.off() # Vero ribo

# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1550, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 15), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.88, "npc"))
# dev.off() # Vero ribo short

# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1400, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 20), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.89, "npc"))
# dev.off() # HBEC rseq

# png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1400, res = 300)
# draw(h, padding = unit(c(2, 5, 8, 19), "mm")) # padding: empty space around heatmap, for legends, row names...
# draw(q_legend, x = unit(0.7, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
# draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.87, "npc"))
# dev.off() # Kyle Exp92 rseq

png(filename = file.path(path, "heatmap", paste0(outputname, ".png")), height = 2000, width = 1200, res = 300)
draw(h, padding = unit(c(2, 5, 8, 16), "mm")) # padding: empty space around heatmap, for legends, row names...
draw(q_legend, x = unit(0.6, "npc"), y = unit(0.95, "npc")) # drawing the legend at (x, y) position
draw(ratio_legend, x = unit(0.9, "npc"), y = unit(0.87, "npc"))
dev.off() # Kyle Exp92 ribo

