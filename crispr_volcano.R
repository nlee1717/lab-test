library(ggplot2)
library(ggthemes)
library(ggrepel)

setwd('C:/Users/NLEE/Documents/Career/Kutulay_Lab/scripts')
prefix = "pre_vs_adpt"
df = read.table('test_pre_vs_adpt.gene_summary.txt', header=T)


# ENRICHMENT ANALYSIS PLOT

g2 <- ggplot(data=df, mapping=aes(x=id, y=-log(pos.score))) +
      geom_point(color="gray50") +
      theme_classic() +
      labs(x='sgRNAs', y='RRA Scores') +
      theme(legend.position = 'none',
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size=6)) +
      geom_point(data=df[1:5,], mapping=aes(x=id, y=-log(pos.score)), color="tomato2", size=1.5) +
      geom_text_repel(data=df[1:5,], aes(id, -log(pos.score), label=id), size=2)
    
png(paste(prefix,'_enrich.png',sep=''), 4, 3, units='in', res=300)
g2
dev.off()






# 
# 
# 
# enrich = data.frame(df$id, -log(df$pos.score))
# # rownames(enrich) = df$id
# colnames(enrich) = c("sgRNAs","RRA_Score")
# 
# #select sgRNAs with high RRA score
# enrich_sig <- enrich[1:5,]
# 
# png("enrich_plot.png", height=400, width=500)
# g <- ggplot(data=enrich, mapping=aes(x=sgRNAs, y=RRA_Score)) + 
# 		geom_point(color="gray40") +
# 		theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold")) +
# 		geom_text(data=enrich_sig, aes(sgRNAs, RRA_Score, label=sgRNAs), vjust = "inward", hjust = "inward") +
# 		# geom_text_repel(data=enrich_sig, aes(sgRNAs, RRA_Score, label=sgRNAs)) +
# 		geom_point(data=enrich_sig, mapping=aes(x=sgRNAs, y=RRA_Score), color="red") 
# 
# # ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
# #         geom_point() + 
# #         theme_minimal() +
# #         geom_text_repel() +
# #         scale_color_manual(values=c("blue", "black", "red")) +
# #         geom_vline(xintercept=c(-0.6, 0.6), col="red") +
# #         geom_hline(yintercept=-log10(0.05), col="red")
# 
# dev.off()
# 
