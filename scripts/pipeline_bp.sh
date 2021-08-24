#!/bin/bash

#
# Generated 08/24/2021 By NaKyung Lee
#
# This script intends to be a data-analysis pipeline for Kutluay lab.
#
#
#
#
#
#
#
# ##############################
# #                            #
# # A. Riboprofiling           #
# #                            #
# ##############################
#
# 1. BBDuk
# 	Input(s): raw fastq file
#	Output(s): trimmed, sorted, separated fastq files by barcodes
#
# 2. Mapping (cellular - STAR / viral - BOWTIE)
#	Input(s): fastq files by barcode
#	Output(s): bam & sam files, log files, count_reduced files (cellular), repaired count files (viral)  
#
# ~ Primary Ananlysis
#
# 3. Read Summary  
#	Input(s): log files
#	Output(s): Reads summary table
#
# 4. Read Length Distribution
#	Input(s): bams, transcriptome bams (cellular)
#	Output(s): RLD tables & plots
#
# 5. Viral Read Plot
#	Input(s): repaired counts
#	Output(s): viral read plot
#
# ~
#
# 6. PCA plot 
#
# 7. tcGSA
#
# 8. edgeR* 
#
# 9. riboTISH
#



#
#
# ##############################
# #                            #
# # B. RNA-sequencing          #
# #                            #
# ##############################
#
# 1. Mapping
#	Input(s): raw fastq files
#	Output(s): bam & sam files, log files, counts
#
# ~ Primary Analysis
#
# 2. Reads Summary Table 
#
# 3. Read Lengths Distribution (viral)
#
# 4. Viral Read Plot
#
# ~
#
# 5. PCA plot (cellular)
#
# 6. tcGSA (cellular)
#
# 7. edgeR
#



#
#
# #############################
# #                           #
# # C. DE genes - edgeR       #               
# #                           #
# #############################
#
# edgeR package is used for pairwise comparison on:
# 	minimum 2 repetition of experiments & 2 sets of conditions
#
# 1. DGE object
# 2. gene ID conversion
# 3. filter & dedupe
# 4. (skipped) PCA
# 5. normalization
# 6. MDS plot
# 7. dispersion estimation
# 8. BCV plot
# 9. pairwise comparison
# 10. Smear Plot
# 11. Volcano Plot
# 12. CPM table
# 13. heatmaps: regular / collapsed / logFC
# 14. heatmap tables
#
# * Notes *
# a. switch plots to ggplot 
# b. unify color palette
#
#

   
