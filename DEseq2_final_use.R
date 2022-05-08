
total24lib <- read.csv('../result_coding2/totalCount_24lib.csv',sep=' ')
head(total24lib)
tail(total24lib)
colnames(total24lib)

# remove lib E2.75_22 lib:
total24lib['E2_75_22_filtered_count.txt'] <- NULL

#filter at least 1 read in at least 1 samples
filter <- apply(total24lib, 1, function(x) length(x[x>0])>=1)
filter
filtered <- total24lib[filter,]
length(filtered)
nrow(filtered)
head(filtered)
tail(row.names(filtered),70)
colnames(total24lib)

tail(filtered)

x1 <- as.factor(c("WT", "Null", "Null", "Null", "WT", "WT", "Het", "WT", "Het", "Null", "WT", "WT", "WT", "WT", "Het", "Het", "WT", "Het", "Het", "Het", "Null", "Null", "Null"))

x2 <- as.factor(c("WT_E2p5", "Null_E2p5", "Null_E2p5", "Null_E2p5","WT_E2p5", "WT_E2p5", "Het_E2p5", "WT_E2p5", "Het_E2p5", "Null_E2p5", "WT_E2p5", "WT_E2p5", "WT_E2p75", "WT_E2p75", "Het_E2p75", "Het_E2p75", "WT_E2p75", "Het_E2p75", "Het_E2p75", "Het_E2p75", "Null_E2p75", "Null_E2p75", "Null_E2p75"))

x1
x2
x <- data.frame(matrix(ncol=2,nrow=23))
x[1] <- x1
x[2] <- x2

x <- x1
library(DESeq2)
colnames(filtered)
gene_filtered <- filtered[grep("ENSM",row.names(filtered)),]
tail(gene_filtered)

coldata <- data.frame(row.names=colnames(gene_filtered), x)
coldata
dds <- DESeqDataSetFromMatrix(countData=gene_filtered, colData=coldata, design=~x)

dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
dds

# Regularized log transformation for clustering/heatmaps, etc. very slow. stop here and go the later steps to do diff analysis.
rld <- rlogTransformation(dds)
ncol(assay(rld))
hist(assay(rld))


library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(x))]
mycols
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
tiff("../result_count2_final/sampleDist_total23lib_merge.tiff", w=1000, h=1000, pointsize=15)
heatmap.2(as.matrix(sampleDists), key=T, density.info="none", trace="none",
          col=colorpanel(100, "blue", "yellow"),
          ColSideColors=mycols[x], RowSideColors=mycols[x],
          margin=c(10, 10), main="Sample Distance Matrix",
          Colv = T, Rowv = T)
dev.off()

#this puts samples as rows, genes as columns 
transpose <- t(assay(rld)) 
transpose_df <- as.data.frame(transpose)
#this is the function that does the PCA
pca.data <- prcomp(transpose_df)
scores = as.data.frame(pca.data$x) 
summary(pca.data)

#2d PCA plot:
pdf("../result_count2_final/RNA23lib_merge.pdf")
DESeq2::plotPCA(rld, intgroup="x")
dev.off()

#=========================================================================
#new PCA plot
#this puts samples as rows, genes as columns 
transpose <- t(assay(rld)) 
transpose_df <- as.data.frame(transpose)
transpose_df
tail(transpose_df)
#this is the function that does the PCA
pca.data <- prcomp(transpose_df)
#pca.data$x
scores = as.data.frame(pca.data$x) 
summary(pca.data) #PC1 and 2: 0.3495  0.09122
PC1 <- 0.1643 
PC2 <- 0.1317

library("ggplot2")
library(dplyr)
library(tibble)
library(ggrepel)
use <- scores[,c(1:2)]
x
use$group <- x
use$stage <- c(rep("E2p5", 12), rep("E2p75", 11))
use$genotype <- c("WT", "Null", "Null", "Null", "WT", "WT", "Het", "WT", "Het", "Null", "WT", "WT", "WT", "WT", "Het", "Het", "WT", "Het", "Het", "Het", "Null", "Null", "Null")

dotSize <- 3
alp <- 0.8
name <- "PCA of 23 libs"
colnames(use) #group, stage, method, poolVSsingle, ps_method
#here can choose which is used to color the plot
#library("ggsci")
ggplot(use, aes(x, y, color=genotype)) +
  geom_point(aes((x = use[,1]) , y = use[,2]), size=dotSize, alpha=alp) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        axis.text.x = element_text(size=8,colour = "black"),
        axis.text.y = element_text(size=8, colour = "black")) + 
  labs(title=paste0(name), x=paste0(colnames(use)[1],": ", PC1), y = paste0(colnames(use)[2], ": ", PC2)) +
  geom_abline(intercept = 0, linetype= "dashed",slope = 0, color = 'gray', size =0.5) +
  geom_vline(linetype= "dashed",xintercept = 0, color = 'gray', size =0.5) +
  coord_fixed()

ggsave(paste0("../result_count2_final/PCA_23lib", "genotype", ".tiff"))
#====================================================================

#=======================================================================
#differential expressed gene lists: not finished
x
#dds
#con1 <- c("Null_E2p5", "Het_E2p5","Null_E2p75", "Het_E2p75","Null_E2p5", "Null_E2p75")
#con2 <- c("WT_E2p5", "WT_E2p5", "WT_E2p75", "WT_E2p75","Het_E2p5","Het_E2p75")
con1 <- c('Null','Het',"Null")
con2 <- c('WT','WT','Het')
length(con1)
for (i in 1:length(con1)) {
  con1[i]
  con2[i]
  res <- results(dds, contrast=c("x", con1[i], con2[i]), alpha=0.1)
  resFilename <- paste0("../result_count2_final/",con1[i], "_vs_", con2[i], "_merge.csv")
  write.csv(res, resFilename)
  
  #below do MA plot
  library(ggplot2)
useForPlot <- data.frame(log10(res$baseMean),res$log2FoldChange,res$padj)
colnames(useForPlot) <- c("log10baseMean", "log2FC", "padj")
useForPlot
up <- useForPlot[which((useForPlot$padj < 0.01) & (useForPlot$log2FC>0)), ]
down <- useForPlot[which((useForPlot$padj < 0.01) & (useForPlot$log2FC<0)), ]
ncol(down)
down
yuplim = 8.0
ydownlim = -8.0

upout <-useForPlot[which(useForPlot$log2FC>= yuplim), ]
if (length(upout$log2FC) >0) {
  upout$log2FC <- yuplim
} 
downout <-useForPlot[which(useForPlot$log2FC<= ydownlim), ]
if (length(downout$log2FC)>0) {
  downout$log2FC <- ydownlim
}


dotSize <- 2
ggplot()+
  geom_point(data = useForPlot, aes(x = useForPlot$log10baseMean, y = useForPlot$log2FC), size = dotSize, alpha=0.1) +
  geom_point(data = up, aes(x = up$log10baseMean, y = up$log2FC), size = dotSize, color="red") + 
  geom_point(data = down, aes(x = down$log10baseMean, y = down$log2FC), size = dotSize, color="blue") +
  labs(title=paste0(con1[i], "_vs_", con2[i]), x="log10baseMean", y = "log2FC") + ylim(-8,8) + geom_point(data = upout, aes(x = upout$log10baseMean, y = upout$log2FC), size = dotSize, color="red", shape=5) + geom_point(data= downout, aes(x = downout$log10baseMean, y = downout$log2FC), size = dotSize, color="blue", shape=5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=0, linetype="dashed", color = "gray", size=1) 

ggsave(paste0("../result_count2_final/",con1[i], "_vs_", con2[i], "_merge.tiff"))

}

