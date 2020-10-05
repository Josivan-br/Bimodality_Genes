library(readr)
library(dplyr)

library(clusterProfiler)
library(org.Hs.eg.db)


setwd("E:/ProjetosGit/bimodality_Genes/Analise_GO/clusterProfie")

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

#
gmpr2_signicative_DESeq2 <- gmpr2_DESeq2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
kctd9_signicative_DESeq2 <- kctd9_DESeq2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
lyrm7_signicative_DESeq2 <- lyrm7_DESeq2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
raver2_signicative_DESeq2 <- raver2_DESeq2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
rrp7bp_signicative_DESeq2 <- rrp7bp_DESeq2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
tm7sf3_signicative_DESeq2 <- tm7sf3_DESeq2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)

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

##
gmpr2_signicative_edgeR <- gmpr2_edgeR %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
kctd9_signicative_edgeR <- kctd9_edgeR %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
lyrm7_signicative_edgeR <- lyrm7_edgeR %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
raver2_signicative_edgeR <- raver2_edgeR %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
rrp7bp_signicative_edgeR <- rrp7bp_edgeR %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
tm7sf3_signicative_edgeR <- tm7sf3_edgeR %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)

######  Genes comuns FD 2

genes_gmpr2 <- gmpr2_signicative_DESeq2[(gmpr2_signicative_DESeq2$genename %in% gmpr2_signicative_edgeR$genename),] 
genes_kctd9 <- kctd9_signicative_DESeq2[(kctd9_signicative_DESeq2$genename %in% kctd9_signicative_edgeR$genename),]
genes_lyrm7 <- lyrm7_signicative_DESeq2[(lyrm7_signicative_DESeq2$genename %in% lyrm7_signicative_edgeR$genename),]
genes_raver2 <-raver2_signicative_DESeq2[(raver2_signicative_DESeq2$genename %in% raver2_signicative_edgeR$genename),]
genes_rrp7bp <-rrp7bp_signicative_DESeq2[(rrp7bp_signicative_DESeq2$genename %in% rrp7bp_signicative_edgeR$genename),]
genes_tm7sf3 <-tm7sf3_signicative_DESeq2[(tm7sf3_signicative_DESeq2$genename %in% tm7sf3_signicative_edgeR$genename),]



#####

genelist_significative <- union(
  genes_gmpr2$genename,  
  union(genes_kctd9$genename, genes_lyrm7$genename)
)
genelist_significative <- union(
  genelist_significative,
  union(genes_raver2$genename,genes_rrp7bp$genename)
)
genelist_significative <- union(
genelist_significative, genes_tm7sf3$genename)

######
significative <- DESeq2 %>% 
  filter(genename %in% genelist_significative)

#

library(pheatmap)

genes_significative2 <- significative

rownames(genes_significative2) <- genes_significative2$genename
genes_significative2$genename <- NULL

genes_significative2 <- as.matrix(genes_significative2)
pheatmap(genes_significative2, scale = "row")


############################## genes FD 2

eggmpr2 = bitr(genes_gmpr2$genename,fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(eggmpr2) <- c("genename","ENTREZID")

egkctd9 = bitr(genes_kctd9$genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(egkctd9) <- c("genename","ENTREZID")

eglyrm7 = bitr(genes_lyrm7$genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(eglyrm7) <- c("genename","ENTREZID")

egraver2 = bitr(genes_raver2$genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(egraver2) <- c("genename","ENTREZID")

egrrp7bp = bitr(genes_rrp7bp$genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(egrrp7bp) <- c("genename","ENTREZID")

egtm7sf3 = bitr(genes_tm7sf3$genename, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(egtm7sf3) <- c("genename","ENTREZID")

listaSamples <- list(GMPR2= eggmpr2$ENTREZID, KCTD9 = egkctd9$ENTREZID, LYRM7= eglyrm7$ENTREZID,
                     RAVER2 = egraver2$ENTREZID, RRP7BP = egrrp7bp$ENTREZID, TM7SF3 = egtm7sf3$ENTREZID)

ck <- compareCluster(geneCluster = listaSamples, fun = "enrichGO", OrgDb="org.Hs.eg.db", ont = "BP")
#ck <- compareCluster(geneCluster = listaSamples, fun = "enrichGO", OrgDb="org.Hs.eg.db", ont = "CC")
#ck <- compareCluster(geneCluster = listaSamples, fun = "enrichGO", OrgDb="org.Hs.eg.db", ont = "MF")

dotplot(ck)

GO <- data.frame(ck)

write.table(GO, "GO_BP_KIRC_enrichGO_FG2_Pv0.01.tsv", quote = F, sep = "\t")
#write.table(GO, "GO_CC_KIRC_enrichGO_FG2_Pv0.01.tsv", quote = F, sep = "\t")
#write.table(GO, "GO_MF_KIRC_enrichGO_FG2_Pv0.01.tsv", quote = F, sep = "\t")

ck <- compareCluster(geneCluster = listaSamples, fun = "enrichKEGG")

dotplot(ck)

VR <- data.frame(ck)
write.table(VR, "VR_KIRC_enrichKEGG_FG3_Pv0.01.csv", quote = F, sep = "\t")

##########################


gmpr2_list <- merge(eggmpr2, genes_gmpr2, by = "genename")
gmpr2_lista <- gmpr2_list[-1]
row.names(gmpr2_lista) <- as.character(gmpr2_lista[,1])
gmpr2_lista <- gmpr2_lista[-1]
gmpr2_lista <- gmpr2_lista[5]

#
kctd9_list <- merge(egkctd9, genes_kctd9, by = "genename")
kctd9_lista <- kctd9_list[-1]
row.names(kctd9_lista) <- as.character(kctd9_lista[,1])
kctd9_lista <- kctd9_lista[-1]
kctd9_lista <- kctd9_lista[5]

#

lyrm7_list <- merge(eglyrm7, genes_lyrm7, by = "genename")
lyrm7_lista <- lyrm7_list[-1]
row.names(lyrm7_lista) <- as.character(lyrm7_lista[,1])
lyrm7_lista <- lyrm7_lista[-1]
lyrm7_lista <- lyrm7_lista[5]

#

raver2_list <- merge(egraver2, genes_raver2, by = "genename")
raver2_lista <- raver2_list[-1]
row.names(raver2_lista) <- as.character(raver2_lista[,1])
raver2_lista <- raver2_lista[-1]
raver2_lista <- raver2_lista[5]

#

rrp7bp_list <- merge(egrrp7bp, genes_rrp7bp, by = "genename")
rrp7bp_lista <- rrp7bp_list[-1]
row.names(rrp7bp_lista) <- as.character(rrp7bp_lista[,1])
rrp7bp_lista <- rrp7bp_lista[-1]
rrp7bp_lista <- rrp7bp_lista[5]

#

tm7sf3_list <- merge(egtm7sf3, genes_tm7sf3, by = "genename")
tm7sf3_lista <- tm7sf3_list[-1]
row.names(tm7sf3_lista) <- as.character(tm7sf3_lista[,1])
tm7sf3_lista <- tm7sf3_lista[-1]
tm7sf3_lista <- tm7sf3_lista[5]

library("pathview")

#Collecting duct acid secretion - hsa04966
# TM7SF3

hsa04966_tm7sf3 <- pathview(gene.data  = tm7sf3_lista,
                            pathway.id = "hsa04966",
                            species    = "hsa",
                            out.suffix = "Collecting duct acid secretion_tm7sf3")


#Vibrio cholerae infection - hsa05110
# TM7SF3

hsa05110_tm7sf3 <- pathview(gene.data  = tm7sf3_lista,
                            pathway.id = "hsa05110",
                            species    = "hsa",
                            out.suffix = "Vibrio cholerae infection_tm7sf3")


## Synaptic vesicle cycle - hsa04721
TM7SF3


hsa04721_tm7sf3 <- pathview(gene.data  = tm7sf3_lista,
                            pathway.id = "hsa04721",
                            species    = "hsa",
                            out.suffix = "Synaptic_Vesicle_Cycle_TM7SF3")


# Mineral absorption - hsa04978
GMPR2


hsa04978_gmpr2 <- pathview(gene.data  = gmpr2_lista,
                            pathway.id = "hsa04978",
                            species    = "hsa",
                            out.suffix = "Mineral absorption_GMPR2")


# Nicotinate and nicotinamide metabolism - hsa00760
GMPR2


hsa00760_gmpr2 <- pathview(gene.data  = gmpr2_lista,
                           pathway.id = "hsa00760",
                           species    = "hsa",
                           out.suffix = "Nicotinate and nicotinamide metabolism_GMPR2")
