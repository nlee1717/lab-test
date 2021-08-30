## Calculates P-site information using riboWaltz and plotting the results. In this case, it uses a one-line SARS-Cov-2 annotation that treats its genome as one transcript. ##

library(riboWaltz)
library(data.table)
library(tools)
library(ggplot2)
library(dplyr)


bam_folder <- "C:/kutluaylab/data/HIV/psite_frame/bam_in_use"
output_folder <- "C:/kutluaylab/data/HIV/psite_frame/MO_outputs/viral"
name <- "RP2_72hpi"

# bam_folder <- "C:/kutluaylab/data/Vero/riboseq/psite_frame/bam_in_use"
# output_folder <- "C:/kutluaylab/data/Vero/riboseq/psite_frame/test_output"
# name <- "Vero_176_6hpi_test"

# custom psite (RiboWaltz function) script #
# source(file = "C:/kutluaylab/scripts/R/psite_frame/psite_edited.R")



# creating annotation #

# annotation_human <- create_annotation(gtfpath = "C:/kutluaylab/data/annotation_files/Homo_sapiens.GRCh37.87_notab.gtf") # human
# annotation_Vero <- create_annotation(gtfpath = "C:/kutluaylab/data/annotation_files/Vero/chlSab2.ncbiRefSeq.gtf") # Vero
# annotation <- create_annotation(gtfpath = "C:/kutluaylab/data/annotation_files/SARS-Cov-2/COVID-19_sepORF1AB_noheader_agat.gtf") # COVID
# annotation <- create_annotation(gtfpath = "C:/kutluaylab/data/annotation_files/SARS-Cov-2/COVID-19_sepORF1AB_noheader_agat_test.gtf") # COVID
# annotation_table <- write.csv(annotation, file = "C:/kutluaylab/data/annotation_files/SARS-Cov-2/COVID-19_riboWaltz_annotation_table.csv")

# save(annotation_human, file = "C:/kutluaylab/data/HIV/psite_frame/annotations/human_annotation.RData")
load(file = "C:/kutluaylab/data/HIV/psite_frame/annotations/human_annotation.RData")


# Custom one-line annotation table #
# annotation_new <- read.csv(file = "C:/kutluaylab/data/annotation_files/SARS-Cov-2/COVID-19_riboWaltz_annotation_table_one_edited.csv",
#                            row.names = 1)
annotation_new <- read.delim(file = "C:/kutluaylab/data/annotation_files/HIV/NL4-3_one_gagpol_edited.txt")
annotation_new <- setDT(annotation_new)


# bam_folder only has one bam in it #
reads_list <- bamtolist(bamfolder = bam_folder, annotation = annotation_new)
# reads_list_human <- bamtolist(bamfolder = "C:/kutluaylab/data/Vero/riboseq/psite_frame/human_bams", annotation = annotation_human)

save(reads_list, file = file.path(output_folder, paste0(name, "_reads_list.RData")))



# read length distribution for getting proper length filter range #
rlen_dist <- rlength_distr(data = reads_list,
                           sample = file_path_sans_ext(dir(bam_folder)))
write.csv(rlen_dist[[1]], file = file.path(output_folder, "RLD", paste0(name, "_rld_table.csv")))

png(filename = file.path(output_folder, "RLD", paste0(name, "_rld.png")),
    width = 1000, height = 500, res = 100)
print(rlen_dist[[2]]) # "Object 'SK393' not found" error
dev.off()



# filtering reads with a length range #
reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_range = 29:34)
# reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "custom", length_filter_vector = 25:33)
# reads_list_filter <- length_filter(data = reads_list, length_filter_mode = "periodicity", periodicity_threshold = 40)

name <- "RP2_72hpi_29_34"
write.csv(reads_list_filter[[1]],
          file = file.path(output_folder, paste0(name, "_reads_filter.csv")))
save(reads_list_filter, file = file.path(output_folder, paste0(name, "_reads_list_filter.RData")))


psite_offset <- psite(reads_list_filter, flanking = 6, extremity = "auto")
# psite_offset <- cusPsite(reads_list_filter)     # cusPsite is custom psite function from psite_edited.R
write.csv(psite_offset,
          file = file.path(output_folder, paste0(name, "_psite_offset.csv")))


reads_psite_list <- psite_info(reads_list_filter, psite_offset)

write.csv(psite_offset, file = file.path(output_folder, paste0(name, "_psite_offset.csv")))



# trying to use cellular offsets for viral reads #

psite_offset_cellular <- read.csv("C:/kutluaylab/data/HIV/psite_frame/MO_outputs/cellular/RP2_cellular_72hpi_29_34_psite_offset.csv",
                                  row.names = 1)
reads_list_filter_table <- read.csv(file.path(output_folder, paste0(name, "_reads_filter.csv")),
                                    row.names = 1)


corrected_offset_5 <- psite_offset_cellular[, c("length", "corrected_offset_from_5")]

# SK293 25-36
# corrected_offset_5$corrected_offset_from_5 <- c(10, 10, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14) # v1
# corrected_offset_5$corrected_offset_from_5 <- c(9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15) # v2
# corrected_offset_5$corrected_offset_from_5 <- c(9, 10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14) # v3
# corrected_offset_5$corrected_offset_from_5 <- c(9, 10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 14) # v4

# SK393 29-34
# corrected_offset_5$corrected_offset_from_5 <- c(11, 12, 12, 13, 13, 13) # v1
# corrected_offset_5$corrected_offset_from_5 <- c(12, 12, 12, 13, 13, 14) # v2
# corrected_offset_5$corrected_offset_from_5 <- c(12, 12, 13, 13, 13, 13) # v3
# corrected_offset_5$corrected_offset_from_5 <- c(12, 13, 13, 13, 13, 13) # v4
# corrected_offset_5$corrected_offset_from_5 <- c(13, 13, 13, 13, 13, 13) # v5

# SK396 29-33
corrected_offset_5 <- corrected_offset_5[1:5, ]
# SK396 29-34
corrected_offset_5$corrected_offset_from_5 <- c(13, 13, 13, 13, 13, 13) # v1
# corrected_offset_5$corrected_offset_from_5 <- c(12, 12, 13, 13, 13, 13) # v2
# corrected_offset_5$corrected_offset_from_5 <- c(12, 13, 13, 13, 13, 13) # v3
# corrected_offset_5$corrected_offset_from_5 <- c(13, 13, 13, 13, 13, 14) # v4

# MO RP2 24hpi
corrected_offset_5 <- corrected_offset_5[1:4, ]



pinfo <- reads_list_filter[[1]]
# pinfo <- reads_list_filter_table


pinfo$psite <- pinfo$end5 +
  corrected_offset_5[match(pinfo$length, corrected_offset_5$length), ]$corrected_offset_from_5

pinfo$psite_from_start <- pinfo$psite - pinfo$cds_start
pinfo$psite_from_stop <- pinfo$psite - pinfo$cds_stop

pinfo$psite_region <- case_when(pinfo$psite_from_start < 0 ~ "5utr",
                                pinfo$psite_from_start >= 0 & pinfo$psite_from_stop <= 0 ~ "cds",
                                pinfo$psite_from_stop > 0 ~ "3utr")

pinfo <- relocate(pinfo, psite, .after = end5)



# pinfo has P-site and other information for each read #
# pinfo <- reads_psite_list[[1]] # the first list element, as there is only one sample
pinfo$frame_from_start <- pinfo$psite_from_start %% 3
pinfo <- pinfo[order(pinfo$psite), ]



gag <- table(pinfo[pinfo$psite >= 336 & pinfo$psite <= 1632, ]$frame_from_start)
pol <- table(pinfo[pinfo$psite >= 1838 & pinfo$psite <= 4640, ]$frame_from_start)


# calculating percentage and writing to file #

gag_perc <- gag["0"] / sum(gag)
pol_perc <- pol["2"] / sum(pol)
gag <- append(gag, gag_perc)
pol <- append(pol, pol_perc)
gag <- format(gag, scientific = F)
pol <- format(pol, scientific = F)
gag[1:3] <- as.integer(gag[1:3])
pol[1:3] <- as.integer(pol[1:3])
names(gag) <- c("0", "1", "2", "gag_perc")
names(pol) <- c("0", "1", "2", "pol_perc")


name <- "RP1_24hpi_2_28_30"

write.csv(data.frame(names(gag), gag),
          file = file.path(output_folder, paste0(name, "_gag.csv")),
          row.names = F)
write.csv(data.frame(names(pol), pol),
          file = file.path(output_folder, paste0(name, "_pol.csv")),
          row.names = F)


write.csv(pinfo, file = file.path(output_folder, paste0(name, "_psite_info.csv")))



name <- "RP2_24hpi_28_34"

pinfo <- read.csv(file.path(output_folder, paste0(name, "_psite_info.csv")), row.names = 1)


# for 176 ORF1AB #
# orf1a <- pinfo[pinfo$psite <= 13467, ]
# orf1b <- pinfo[pinfo$psite >= 13468 & pinfo$psite <= 21555, ]
# orf1a_table <- data.frame(table(orf1a$frame_from_start))
# orf1b_table <- data.frame(table(orf1b$frame_from_start))
#
# colnames(orf1a_table) <- c("frame", "count")
# colnames(orf1b_table) <- c("frame", "count")
#
# write.csv(orf1a_table, file = file.path(output_folder, paste0(name, "_whole_1a.csv")))
# write.csv(orf1b_table, file = file.path(output_folder, paste0(name, "_whole_1b.csv")))




## Attempt to plot P-site counts and their frames across the entire ORF1AB region ##

# orf1ab <- pinfo[pinfo$psite <= 21555, ]
# orf1ab <- orf1ab[order(orf1ab$psite), ]

# orf1ab_table <- data.frame(table(orf1ab$psite))
# colnames(orf1ab_table) <- c("psite", "count")
# orf1ab_table$psite <- as.integer(as.vector(orf1ab_table$psite))
# orf1ab_table$frame <- as.integer((orf1ab_table$psite + 1) %% 3)
#
#
# p1 <- (ggplot(data = orf1ab_table, mapping = aes(x = psite, y = count, fill = as.factor(frame)))
#   + geom_col(width = 1)
#
#   + scale_x_continuous(name = "P-site Location"
#                        # breaks = c(13450, 13468, 13500),
#                        # labels = c("13450", "Frameshifting\nSite", "13500")
#                       )
#   + scale_y_continuous(name = "Count")
#
#   + scale_fill_discrete(name = "Frame")
#
#   + theme_classic()
#   + theme(legend.position = c(0.85, 0.9),
#           legend.title = element_text(size = 5),
#           axis.title = element_text(size = 5))
#
#   + guides(fill = guide_legend(nrow = 1))
# )
#
# png(filename = file.path(output_folder, paste0(name, "_test_plot.png")), width = 4000, height = 1500, res = 1000)
# print(p1)
# dev.off()
#
# png(filename = file.path(output_folder, paste0(name, "_orf1ab_plot.png")), width = 4000, height = 1500, res = 500)
# print(p1)
# dev.off()





## Selecting the frameshift region ##

fs <- pinfo[pinfo$psite >= 1500 & pinfo$psite < 2000, ]
# fs <- fs[order(fs$psite), ]
write.csv(fs, file = file.path(output_folder, paste0(name, "_frameshift.csv")))

name <- "176_24hpi_collapsed"
fs <- read.csv(file.path(output_folder, paste0(name, "_frameshift.csv")), row.names = 1)



## P-site frame distribution of orf1a and orf1b regions ##

fs_gag <- fs[fs$psite <= 1629, ]
fs_pol <- fs[fs$psite >= 1631, ]
fs1a_table <- data.frame(table(fs1a$frame_from_start))
fs1b_table <- data.frame(table(fs1b$frame_from_start))
colnames(fs1a_table) <- c("frame", "count")
colnames(fs1b_table) <- c("frame", "count")

write.csv(fs1a_table, file = file.path(output_folder, paste0(name, "_frameshift_1a.csv")))
write.csv(fs1b_table, file = file.path(output_folder, paste0(name, "_frameshift_1b.csv")))




## PLOTTING ##
## plots count of all p-sites across a genome region as height of columns,
## coloring the columns based on their frames ##

fs_short <- fs[, c("psite", "frame_from_start")]
fs_short$frame_from_start <- as.integer(fs_short$frame_from_start)


# breaks <- sort(c(pretty(fs_short$psite), 1631, 1838))
# labels <- sub("1631", "Frameshifting\nsite", as.character(breaks))
# labels <- sub("1838", "Gag stop\ncodon", labels)

breaks <- c(1500, 1631, 1700, 1838, 1900)
labels <- c("1500", "Frameshifting\nsite", "1700", "Gag stop\ncodon", "1900")


p <- (ggplot(data = fs_short, mapping = aes(x = psite, fill = as.factor(frame_from_start)))
      + geom_histogram(binwidth = 1) # can also use geom_col when we have count of each p-site

      + scale_x_continuous(name = "P-site Location",
                           breaks = breaks,
                           labels = labels)
      + scale_y_continuous(name = "Count")

      + scale_fill_discrete(name = "Frame")

      + theme_classic()
      + theme(legend.position = c(0.85, 0.9),
              axis.title = element_text(size = 13))

      + guides(fill = guide_legend(nrow = 1))
)

png(filename = file.path(output_folder, paste0(name, "_plot.png")), width = 2000, height = 1500, res = 300)
print(p)
dev.off()
