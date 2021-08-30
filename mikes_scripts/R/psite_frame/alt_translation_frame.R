library(riboWaltz)
library(data.table)


bam_folder <- "C:/kutluaylab/data/Vero/riboseq/psite_frame/bam_in_use"
output_folder <- "C:/kutluaylab/data/Vero/riboseq/psite_frame/RP9_outputs_2"
name <- "RP9_6hpi_HHT_collapsed"

annotation_new <- read.delim(file = "C:/kutluaylab/data/Vero/riboseq/psite_frame/annotation/COVID-19_riboWaltz_annotation_table_regions.csv",
                            row.names = 1)
annotation_new <- setDT(annotation_new)

reads_list <- bamtolist(bamfolder = bam_folder, annotation = annotation_new)

reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 25:33)

psite_offset <- psite(reads_list_filter, flanking = 6, extremity = "auto")

reads_psite_list <- psite_info(reads_list_filter, psite_offset)

write.csv(psite_offset, file = file.path(output_folder, paste0(name, "_psite_offset.csv")))


pinfo <- reads_psite_list[[1]]
pinfo$frame_from_start <- pinfo$psite_from_start %% 3

write.csv(pinfo, file = file.path(output_folder, paste0(name, "_psite_info.csv")))



folder <- "C:/kutluaylab/data/Vero/riboseq/psite_frame/RP9_outputs"
name <- "RP9_24hpi_HHT_collapsed"

pinfo <- read.csv(file.path(folder, paste0(name, "_psite_info.csv")), row.names = 1)
psite <- pinfo$psite



# numbers are converted from 0-based coordinates from ribo-TISH table to 1-based coordinates

# S(6) #
# alt_trans <- "s6"
# region <- 21868:21924
# anno_tis <- 21563

# E(6) #
# alt_trans <- "e6"
# region <- 26389:26472
# anno_tis <- 26245

# M(6) #
# alt_trans <- "m6"
# region <- 26711:26872
# anno_tis <- 26523

# M(24) #
alt_trans <- "m24"
region <- 26658:27191
anno_tis <- 26523


alt_frame_num <- (region[1] - anno_tis) %% 3 # the frame of alternative translation site, counting from the annotated TIS

psite_region <- psite[psite %in% region]
frame <- (psite_region - anno_tis) %% 3
frame_bias <- table(frame)
write.csv(data.frame(frame_bias), file = file.path(folder, paste(name, alt_trans, "alt_trans_table.csv",
                                                                 sep = "_")))



upstream_region <- (region[1] - 99):(region[1] - 1)

psite_upstream <- psite[psite %in% upstream_region]
frame_upstream <- (psite_upstream - anno_tis) %% 3
frame_bias_upstream <- table(frame_upstream)




## Plotting ##

psite_frame_table <- data.frame(c(psite_upstream, psite_region), as.integer(c(frame_upstream, frame)))
colnames(psite_frame_table) <- c("psite", "frame")

write.csv(psite_frame_table, file = file.path(folder, paste(name, alt_trans, "alt_trans_psite_frame.csv",
                                                            sep = "_")))



## Reading the psite-frame table to plot again ##

folder <- "C:/kutluaylab/data/Vero/riboseq/psite_frame/RP9_outputs"
name <- "RP9_24hpi_HHT_collapsed"
alt_trans <- "m24"
region <- 26658:27191

psite_frame_table <- read.csv(file = file.path(folder, paste(name, alt_trans, "alt_trans_psite_frame.csv",
                                                             sep = "_")),
                              row.names = 1)


## Custom x-axis breaks and labels ##
breaks <- sort(c(pretty(psite_frame_table$psite), region[1]))
labels <- sub(as.character(region[1]), "Alt.\nTIS", as.character(breaks))


p <- (ggplot(data = psite_frame_table, mapping = aes(x = psite, fill = as.factor(frame)))
  + geom_histogram(binwidth = 1)

  + scale_x_continuous(name = "P-site Location",
                       breaks = breaks,
                       labels = labels
                      )
  + scale_y_continuous(name = "Count")

  + scale_fill_discrete(name = "Frame")

  + theme_classic()
  + theme(legend.position = c(0.2, 0.9),
          axis.title = element_text(size = 13))

  + guides(fill = guide_legend(nrow = 1))
)



png(filename = file.path(folder, paste(name, alt_trans, "alt_trans_plot_withTIS.png",
                                       sep = "_")),
    width = 2000, height = 1500, res = 300)
print(p)
dev.off()

