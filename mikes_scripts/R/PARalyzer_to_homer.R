## attempt to convert PARalyzer cluster output to HOMER peak format.

par <- read.csv(file = "C:/kutluaylab/data/CLIP_test/Exp492_cluster.csv")


homer <- data.frame(par$ClusterID, par$Chromosome, par$ClusterStart, par$ClusterEnd, par$Strand)


# par_tab <- NULL
#
# for (i in seq(ncol(par))) {
#   par_tab <- cbind(par_tab, par[[i]])
# }

write.table(homer,
            file = "C:/kutluaylab/data/CLIP_test/Exp492_cluster_homer.txt",
            quote = F,
            sep = "\t",
            col.names = F,
            row.names = F)

