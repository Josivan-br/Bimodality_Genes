library(readr)
library(dplyr)

library(clusterProfiler)
library(org.Hs.eg.db)


setwd("E:/ProjetosGit/bimodality_Genes/Analise_s-score")

#################### PRIMEIRO FILTRO FD 2 ############################

####  DESeq2 FD2

columns <- c("genename", "baseMean", "baseMeanA",	"baseMeanB", "foldChange",	"log2FoldChange",	"pval",	"FDR")

gmpr2_DESeq2 <- read_tsv(file = "./gmpr2_DESeq2.tsv", col_names = columns, skip = 1)
kctd9_DESeq2 <- read_tsv(file = "./kctd9_DESeq2.tsv", col_names = columns, skip = 1)
lyrm7_DESeq2 <- read_tsv(file = "./lyrm7_DESeq2.tsv", col_names = columns, skip = 1)
raver2_DESeq2 <- read_tsv(file = "./raver2_DESeq2.tsv", col_names = columns, skip = 1)
rrp7bp_DESeq2 <- read_tsv(file = "./rrp7bp_DESeq2.tsv", col_names = columns, skip = 1)
tm7sf3_DESeq2 <- read_tsv(file = "./tm7sf3_DESeq2.tsv", col_names = columns, skip = 1)

#

DESeq2 <- merge(gmpr2_DESeq2[, c("genename", "log2FoldChange")], 
                kctd9_DESeq2[, c("genename", "log2FoldChange")], 
                by = "genename")

DESeq2 <- merge(DESeq2, lyrm7_DESeq2[, c("genename", "log2FoldChange")], by = "genename")
DESeq2 <- merge(DESeq2, raver2_DESeq2[, c("genename", "log2FoldChange")], by = "genename")
DESeq2 <- merge(DESeq2, rrp7bp_DESeq2[, c("genename", "log2FoldChange")], by = "genename")
DESeq2 <- merge(DESeq2, tm7sf3_DESeq2[, c("genename", "log2FoldChange")], by = "genename")

colnames(DESeq2) <- c("genename", "GMPR2", "KCTD9", "LYRM7", "RAVER2", "RRP7BP", "TM7SF3")


###### edgeR  FD2

columns <- c("genename", "baseMean", "baseMeanA",	"baseMeanB", "foldChange",	"log2FoldChange",	"pval",	"FDR")

gmpr2_edgeR <- read_tsv(file = "./gmpr2_edgeR.tsv", col_names = columns, skip = 1)
kctd9_edgeR <- read_tsv(file = "./kctd9_edgeR.tsv", col_names = columns, skip = 1)
lyrm7_edgeR <- read_tsv(file = "./lyrm7_edgeR.tsv", col_names = columns, skip = 1)
raver2_edgeR <- read_tsv(file = "./raver2_edgeR.tsv", col_names = columns, skip = 1)
rrp7bp_edgeR <- read_tsv(file = "./rrp7bp_edgeR.tsv", col_names = columns, skip = 1)
tm7sf3_edgeR <- read_tsv(file = "./tm7sf3_edgeR.tsv", col_names = columns, skip = 1)

#

edgeR <- merge(gmpr2_edgeR[, c("genename", "log2FoldChange")], 
               kctd9_edgeR[, c("genename", "log2FoldChange")], 
               by = "genename")

edgeR <- merge(edgeR, lyrm7_edgeR[, c("genename", "log2FoldChange")], by = "genename")
edgeR <- merge(edgeR, raver2_edgeR[, c("genename", "log2FoldChange")], by = "genename")
edgeR <- merge(edgeR, rrp7bp_edgeR[, c("genename", "log2FoldChange")], by = "genename")
edgeR <- merge(edgeR, tm7sf3_edgeR[, c("genename", "log2FoldChange")], by = "genename")

colnames(edgeR) <- c("genename", "GMPR2", "KCTD9", "LYRM7", "RAVER2", "RRP7BP", "TM7SF3")

######  Genes comuns FD 2

genes_gmpr2 <- gmpr2_edgeR[(gmpr2_edgeR$genename %in% gmpr2_DESeq2$genename),] 
genes_kctd9 <- kctd9_edgeR[(kctd9_edgeR$genename %in% kctd9_DESeq2$genename),]
genes_lyrm7 <- lyrm7_edgeR[(lyrm7_edgeR$genename %in% lyrm7_DESeq2$genename),]
genes_raver2 <-raver2_edgeR[(raver2_edgeR$genename %in% raver2_DESeq2$genename),]
genes_rrp7bp <-rrp7bp_edgeR[(rrp7bp_edgeR$genename %in% rrp7bp_DESeq2$genename),]
genes_tm7sf3 <-tm7sf3_edgeR[(tm7sf3_edgeR$genename %in% tm7sf3_DESeq2$genename),]

#### S-SCORE ######

s_score <- read.table("genestats.KIRC.all.txt", sep = "\t", header = T, stringsAsFactors = F)

score_gmpr2 <- s_score[(s_score$genename %in% genes_gmpr2$genename),]
score_kctd9 <- s_score[(s_score$genename %in% genes_kctd9$genename),]
score_lyrm7 <- s_score[(s_score$genename %in% genes_lyrm7$genename),]
score_raver2 <- s_score[(s_score$genename %in% genes_raver2$genename),]
score_rrp7bp <- s_score[(s_score$genename %in% genes_rrp7bp$genename),]
score_tm7sf3 <- s_score[(s_score$genename %in% genes_tm7sf3$genename),]

###
library(dplyr)

score_gmpr2_signicative <- score_gmpr2 %>% filter(abs(score) >= 3)
score_kctd9_signicative <- score_kctd9 %>% filter(abs(score) >= 3)
score_lyrm7_signicative <- score_lyrm7 %>% filter(abs(score) >= 3)
score_raver2_signicative <- score_raver2 %>% filter(abs(score) >= 3)
score_rrp7bp_signicative <- score_rrp7bp %>% filter(abs(score) >= 3)
score_tm7sf3_signicative <- score_tm7sf3 %>% filter(abs(score) >= 3)


genelist_significative <- union(
  score_gmpr2_signicative$genename, 
  union(score_kctd9_signicative$genename, score_lyrm7_signicative$genename)
)
genelist_significative <- union(
  genelist_significative,
  union(score_raver2_signicative$genename,score_tm7sf3_signicative$genename)
)
genelist_significative <- union(
  genelist_significative, score_rrp7bp_signicative$genename)

###

Genes_Score <- rbind(score_gmpr2_signicative[,c("genename", "score")], 
                     score_kctd9_signicative[,c("genename","score")])

Genes_Score <- rbind(Genes_Score, score_lyrm7_signicative[, c("genename", "score")])
Genes_Score <- rbind(Genes_Score, score_raver2_signicative[, c("genename", "score")])
Genes_Score <- rbind(Genes_Score, score_rrp7bp_signicative[, c("genename", "score")])
Genes_Score <- rbind(Genes_Score, score_tm7sf3_signicative[, c("genename", "score")])
Genes_Score <- Genes_Score[!(duplicated(Genes_Score$genename)),]

Genes_Score <- merge(DESeq2,Genes_Score, by="genename")
Genes_Score <- Genes_Score[(complete.cases(Genes_Score)),]

write.table(Genes_Score, "Genes_Score_Pv0.01.tsv", quote = F, sep = "\t")

#

library(pheatmap)

genes_significative2 <- Genes_Score

rownames(genes_significative2) <- genes_significative2$genename
genes_significative2$genename <- NULL

genes_significative2 <- as.matrix(genes_significative2)
pheatmap(genes_significative2, scale = "row")

