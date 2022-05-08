# analyze ChIP-seq consensus peaks.
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html
# install done
# BiocManager::install("ChIPseeker")
# BiocManager::install("EnsDb.Mmusculus.v79")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")

# Load libraries
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#Team BC, Maintainer BP (2019). TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package for TxDb object(s). R package version 3.4.7.
library(EnsDb.Mmusculus.v79)
#Rainer J (2017). EnsDb.Mmusculus.v79: Ensembl based annotation package. R package version 2.99.0.
library(clusterProfiler)
citation("clusterProfiler")
# Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu x, Liu S, Bo X, Yu G (2021). “clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.” The Innovation, 2(3), 100141. doi: 10.1016/j.xinn.2021.100141.

#Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), 284-287. doi: 10.1089/omi.2011.0118.
library(AnnotationDbi)
# Pagès H, Carlson M, Falcon S, Li N (2021). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.56.2, https://bioconductor.org/packages/AnnotationDbi.
library(org.Mm.eg.db)
# Carlson M (2019). org.Mm.eg.db: Genome wide annotation for Mouse. R package version 3.8.2.

library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%


# Load data
# bed files must have chr.
samplefiles <- list.files("/Users/Di/Documents/J_lab/dis3zygoticNew/consensusPeaksChIP", pattern= "chr.bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Pou6f1")

samplefiles

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

peakAnnoList

plotAnnoBar(peakAnnoList)

#Distribution of TF-binding loci relative to TSS
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")


pou6f1_annot <- data.frame(peakAnnoList[["Pou6f1"]]@anno)
pou6f1_annot

# Get the entrez IDs
entrez <- pou6f1_annot$geneId
entrez
keytypes(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Mmusculus.v79,
                                         keys = entrez,
                                         columns = c("GENENAME",'GENEID'),
                                         keytype = "ENTREZID")

# Change IDs to character type to merge
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

annotations_edb
# Write to file
pou6f1_annot

pou6f1_annot %>% 
  left_join(annotations_edb, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="../result_ChIP-seq/pou6f1_anno_peaks.csv", sep="\t", quote=F, row.names=F)


# Run GO enrichment analysis 
ego <- enrichGO(gene = entrez, 
                keyType = "ENTREZID", 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "../result_ChIP-seq/go_analysis.csv")
cluster_summary

# Dotplot visualization
dotplot(ego, showCategory=50)


# do Go enrichment of filtered genes: only choose promoter associated genes:
entrez2 <- pou6f1_annot[pou6f1_annot['annotation']=='Promoter',]$geneId
entrez2

ego2 <- enrichGO(gene = entrez2, 
                keyType = "ENTREZID", 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego2)
write.csv(cluster_summary, "../result_ChIP-seq/go_analysisOfPromoterOccupiedOnly.csv")
cluster_summary


