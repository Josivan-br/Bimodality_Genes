library(readr)
library(dplyr)

library(clusterProfiler)
library(org.Hs.eg.db)


setwd("E:/AmostrasKIRC/KIRC_CONT/clusterProfie")

#load("E:/AmostrasKIRC/KIRC_CONT/clusterProfie/ViasRegulasionKIRC.RData")

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
gmpr2_signicative_DESeq2 <- gmpr2_DESeq2 %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
kctd9_signicative_DESeq2 <- kctd9_DESeq2 %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
lyrm7_signicative_DESeq2 <- lyrm7_DESeq2 %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
raver2_signicative_DESeq2 <- raver2_DESeq2 %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
rrp7bp_signicative_DESeq2 <- rrp7bp_DESeq2 %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
tm7sf3_signicative_DESeq2 <- tm7sf3_DESeq2 %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)

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
gmpr2_signicative_edgeR <- gmpr2_edgeR %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
kctd9_signicative_edgeR <- kctd9_edgeR %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
lyrm7_signicative_edgeR <- lyrm7_edgeR %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
raver2_signicative_edgeR <- raver2_edgeR %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
rrp7bp_signicative_edgeR <- rrp7bp_edgeR %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)
tm7sf3_signicative_edgeR <- tm7sf3_edgeR %>% filter(abs(log2FoldChange) >= 2, pval <= 0.01)

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

#########  primeiro método de amostras 

genes_gmpr2 <- genes_gmpr2[c(1,6)] 
genes_kctd9 <- genes_kctd9[c(1,6)]
genes_lyrm7 <- genes_lyrm7[c(1,6)]
genes_raver2 <- genes_raver2[c(1,6)]
genes_rrp7bp <- genes_rrp7bp[c(1,6)]
genes_tm7sf3 <- genes_tm7sf3[c(1,6)]

##

genes_gmpr2$gene <- "GMPR2" 
genes_kctd9$gene <- "KCDT9"
genes_lyrm7$gene <- "LYRM7"
genes_raver2$gene <- "RAVER2"
genes_rrp7bp$gene <- "RRP7BP"
genes_tm7sf3$gene <- "TM7SF3"



########  segundo método de amostras

significative2 <- significative

genes_gmpr2 <- significative2[c(1,2)] 
colnames(genes_gmpr2) <-c("genename","log2FoldChange")

genes_kctd9 <- significative2[c(1,3)]
colnames(genes_kctd9) <-c("genename","log2FoldChange")

genes_lyrm7 <- significative2[c(1,4)]
colnames(genes_lyrm7) <-c("genename","log2FoldChange")

genes_raver2 <- significative2[c(1,5)]
colnames(genes_raver2) <-c("genename","log2FoldChange")

genes_rrp7bp <- significative2[c(1,6)]
colnames(genes_rrp7bp) <-c("genename","log2FoldChange")

genes_tm7sf3 <- significative2[c(1,7)]
colnames(genes_tm7sf3) <-c("genename","log2FoldChange")

##

genes_gmpr2$gene <- "GMPR2" 
genes_kctd9$gene <- "KCDT9"
genes_lyrm7$gene <- "LYRM7"
genes_raver2$gene <- "RAVER2"
genes_rrp7bp$gene <- "RRP7BP"
genes_tm7sf3$gene <- "TM7SF3"



##
significative2 <- rbind(genes_gmpr2,genes_kctd9,genes_lyrm7,
                     genes_raver2,genes_rrp7bp,genes_tm7sf3)


################  ##########################

library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

if(!require('UpSetR')) { install.packages('UpSetR') }

### Comando para criar um data.frame
genes <-  c("GMPR2", "KCDT9", "LYRM7", "RAVER2", "RRP7BP", "TM7SF3")

up_set <- data.frame(matrix(nrow=length(unique(significative2$genename)), ncol= 6, dimnames = list(unique(significative2$genename), genes) )) # cria um data.frame vazio
up_set[, genes] <- 0


for(g in colnames(up_set)){
  up_set[significative2$genename[significative2$gene == g], g] <-  significative2$log2FoldChange[significative2$gene == g]
}

colnames(up_set) <- toupper(colnames(up_set))

# Plot regulon activity profiles
pdf(file="Heatmap_pheatmap.pdf",width = 11, height = 8)
print(
  pheatmap(up_set, 
           main="Bimodal Genes",
           #annotation_col = hypoxia_regact, 
           show_colnames = TRUE, #annotation_legend = TRUE, 
           clustering_method = "ward.D2", fontsize_row = 5,
           clustering_distance_rows = "correlation")
)
dev.off()


# Another way to heatmap plot

df_molten_dat <- melt(as.matrix(up_set))

names(df_molten_dat)[c(1:2)] <- c("Genes", "Gene_peak")

df_molten_dat$value[df_molten_dat$value == 0] <- NA

# Horizontal
pdf(file="Heatmap_genes_horizontal.pdf",width = 11, height = 8)
print(
  ggplot(data = df_molten_dat,
         aes(x = Genes, y = Gene_peak, fill = value)) + 
    geom_raster() +
    xlab("Genes (logFoldChange)") +
    ylab("Bimodal Genes") + 
    #scale_fill_brewer(palette="RdYlBu", na.value="yellow") +
    scale_fill_gradientn(colours = brewer.pal(7, "RdYlBu"),  na.value="yellow") +
    scale_colour_manual(values=NA) + 
    theme(axis.text.x =  element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(hjust = 1)) +
    labs(fill = "log2FC", title = "Heatmap of Bimodal Genes") + 
    geom_point(data = df_molten_dat, aes(size="", shape = NA), colour = "yeloow") +
    guides(size=guide_legend("NA", override.aes=list(shape=15, size = 7)))
)
dev.off()


# Vertical
pdf(file="Heatmap_genes_vertical.pdf",width = 11, height = 8)
print(
  ggplot(data = df_molten_dat,
         aes(x = Gene_peak, y = Genes, fill = value)) + 
    geom_raster() +
    xlab("Bimodal Genes") +
    ylab("Genes (logFoldChange)") + 
    #scale_fill_brewer(palette="RdYlBu", na.value="yellow") +
    scale_fill_gradientn(colours = brewer.pal(7, "RdYlBu"),  na.value="yellow") +
    scale_colour_manual(values=NA) + 
    theme(axis.text.y =  element_text(hjust = 1, size = 5),
          axis.text.x = element_text(hjust = 1)) +
    labs(fill = "log2FC", title = "Heatmap of Bimodal Genes") +
    geom_point(data = df_molten_dat, aes(size="", shape = NA), colour = "yellow") +
    guides(size=guide_legend("NA", override.aes=list(shape=15, size = 7)))
)
dev.off()




############### ########################

genes_gmpr2 <- genes_gmpr2[c(1,6)] 
genes_kctd9 <- genes_kctd9[c(1,6)]
genes_lyrm7 <- genes_lyrm7[c(1,6)]
genes_raver2 <- genes_raver2[c(1,6)]
genes_rrp7bp <- genes_rrp7bp[c(1,6)]
genes_tm7sf3 <- genes_tm7sf3[c(1,6)]

##

genes_gmpr2$gene <- "GMPR2" 
genes_kctd9$gene <- "KCDT9"
genes_lyrm7$gene <- "LYRM7"
genes_raver2$gene <- "RAVER2"
genes_rrp7bp$gene <- "RRP7BP"
genes_tm7sf3$gene <- "TM7SF3"



##
significative <- rbind(genes_gmpr2,genes_kctd9,genes_lyrm7,
                     genes_raver2,genes_rrp7bp,genes_tm7sf3)



if(!require('UpSetR')) { install.packages('UpSetR') }
genes <- unique(significative$gene)


### Comando para criar um data.frame

up_set <- data.frame(matrix(nrow=length(unique(significative$genename)), ncol= 6, dimnames = list(unique(significative$genename), genes) )) # cria um data.frame vazio
up_set[, genes] <- 0

### Preenchimento do data.frame
for(g in colnames((up_set))){
  up_set[significative$genename[significative$gene == g], g] <-  significative$log2FoldChange[significative$gene == g]
}





#
genes_significative <- DESeq2 %>% 
  filter(genename %in% genelist_significative)

rownames(genes_significative) <- genes_significative$genename
genes_significative$genename <- NULL

genes_significative <- as.matrix(genes_significative)
pheatmap(genes_significative, scale = "row")

write.table(genes_significative, "genes_significative_FD2.tsv", quote = F, sep = "\t")

a <- genes_significative$genename
write.table(a, "list_genes_significative_FD2.tsv", quote = F, sep = "\t")


### Selection genes >=3

signicative_gmpr2 <- genes_gmpr2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
signicative_kctd9 <- genes_kctd9 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
signicative_lyrm7 <- genes_lyrm7 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
signicative_raver2 <- genes_raver2 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
signicative_rrp7bp <- genes_rrp7bp %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)
signicative_tm7sf3 <- genes_tm7sf3 %>% filter(abs(log2FoldChange) >= 3, pval <= 0.01)

##

signicative_gmpr2 <- signicative_gmpr2[c(1,6)] 
signicative_kctd9 <- signicative_kctd9[c(1,6)]
signicative_lyrm7 <- signicative_lyrm7[c(1,6)]
signicative_raver2 <- signicative_raver2[c(1,6)]
#signicative_rrp7bp <- signicative_rrp7bp[c(1,6)]
signicative_tm7sf3 <- signicative_tm7sf3[c(1,6)]

##

signicative_gmpr2$gene <- "GMPR2" 
signicative_kctd9$gene <- "KCTD9"
signicative_lyrm7$gene <- "LYRM7"
signicative_raver2$gene <- "RAVER2"
#signicative_rrp7bp$gene <- "RRP7BP"
signicative_tm7sf3$gene <- "TM7SF3"



##
signicative2 <- rbind(signicative_gmpr2,signicative_kctd9,signicative_lyrm7,
                     signicative_raver2,signicative_tm7sf3)


################  ##########################

library(pheatmap)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


if(!require('UpSetR')) { install.packages('UpSetR') }
genes <- unique(signicative$gene)


### Comando para criar um data.frame

up_set <- data.frame(matrix(nrow=length(unique(signicative$genename)), ncol= 5, dimnames = list(unique(signicative$genename), genes) )) # cria um data.frame vazio
up_set[, genes] <- 0

### Preenchimento do data.frame
for(g in colnames((up_set))){
  up_set[signicative$genename[signicative$gene == g], g] <-  signicative$log2FoldChange[signicative$gene == g]
}

colnames(up_set) <- toupper(colnames(up_set))

# Plot regulon activity profiles
pdf(file="Heatmap_pheatmap.pdf",width = 11, height = 8)
print(
  pheatmap(up_set, 
           main="Bimodal Genes",
           #annotation_col = hypoxia_regact, 
           show_colnames = TRUE, #annotation_legend = TRUE, 
           clustering_method = "ward.D2", fontsize_row = 5,
           clustering_distance_rows = "correlation")
)
dev.off()


# Another way to heatmap plot

df_molten_dat <- melt(as.matrix(up_set))

names(df_molten_dat)[c(1:2)] <- c("Genes", "Gene_peak")

df_molten_dat$value[df_molten_dat$value == 0] <- NA


# Horizontal
pdf(file="Heatmap_genes_horizontal.pdf",width = 11, height = 8)
print(
  ggplot(data = df_molten_dat,
         aes(x = Genes, y = Gene_peak, fill = value)) + 
    geom_raster() +
    xlab("Genes (logFoldChange)") +
    ylab("Bimodal Genes") + 
    #scale_fill_brewer(palette="RdYlBu", na.value="yellow") +
    scale_fill_gradientn(colours = brewer.pal(7, "RdYlBu"),  na.value="yellow") +
    scale_colour_manual(values=NA) + 
    theme(axis.text.x =  element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(hjust = 1)) +
    labs(fill = "log2FC", title = "Heatmap of Bimodal Genes") + 
    geom_point(data = df_molten_dat, aes(size="", shape = NA), colour = "yelow") +
    guides(size=guide_legend("NA", override.aes=list(shape=15, size = 7)))
)
dev.off()


# Vertical
pdf(file="Heatmap_genes_vertical.pdf",width = 11, height = 8)
print(
  ggplot(data = df_molten_dat,
         aes(x = Gene_peak, y = Genes, fill = value)) + 
    geom_raster() +
    xlab("Bimodal Genes") +
    ylab("Genes (logFoldChange)") + 
    #scale_fill_brewer(palette="RdYlBu", na.value="yellow") +
    scale_fill_gradientn(colours = brewer.pal(7, "RdYlBu"),  na.value="yellow") +
    scale_colour_manual(values=NA) + 
    theme(axis.text.y =  element_text(hjust = 1, size = 5),
          axis.text.x = element_text(hjust = 1)) +
    labs(fill = "log2FC", title = "Heatmap of Bimodal Genes") +
    geom_point(data = df_molten_dat, aes(size="", shape = NA), colour = "yellow") +
    guides(size=guide_legend("NA", override.aes=list(shape=15, size = 7)))
)
dev.off()




############### ########################


significative <- union(
  signicative_gmpr2$genename, 
  union(signicative_kctd9$genename, signicative_lyrm7$genename)
)
significative <- union(
  significative,
  union(signicative_raver2$genename,signicative_rrp7bp$genename)
)
significative <- union(
  significative, signicative_tm7sf3$genename)

#



significative <- genes_significative %>% 
  filter(genename %in% significative)


write.table(significative, "genes_significative_FD4.tsv", quote = F, sep = "\t")

library(pheatmap)

significative2 <- significative

rownames(significative2) <- significative2$genename
significative2$genename <- NULL

significative2 <- as.matrix(significative2)
pheatmap(significative2, scale = "row")





############################## Processos GO -  genes FD 2

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
write.table(VR, "VR_KIRC_enrichKEGG_FG2_Pv0.01.csv", quote = F, sep = "\t")

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


## RAS signaling pathway - hsa04014
# RRP7BP

hsa04014_rrp7bp <- pathview(gene.data  = rrp7bp_lista,
                            pathway.id = "hsa04014",
                            species    = "hsa",
                            out.suffix = "RAS_rrp7bp_signaling")


## Neuroactive ligand-receptor interaction - hsa04080
# KCTD9 e TM7SF3

hsa04080_kctd9 <- pathview(gene.data  = kctd9_lista,
                           pathway.id = "hsa04080",
                           species    = "hsa",
                           out.suffix = "kctd9_Neuroactive_ligand-receptor")

hsa04080_tm7sf3 <- pathview(gene.data  = tm7sf3_lista,
                            pathway.id = "hsa04080",
                            species    = "hsa",
                            out.suffix = "tm7sf3_Neuroactive_ligand-receptor")




## Collecting duct acid secretion - hsa04966
# KCTD9, RAVER2, TM7SF3

hsa04966_kctd9 <- pathview(gene.data  = kctd9_lista,
                            pathway.id = "hsa04966",
                            species    = "hsa",
                            out.suffix = "kctd9_Collecting_duct_acid_secretion")

hsa04966_raver2 <- pathview(gene.data  = raver2_lista,
                           pathway.id = "hsa04966",
                           species    = "hsa",
                           out.suffix = "raver2_Collecting_duct_acid_secretion")

hsa04966_tm7sf3 <- pathview(gene.data  = tm7sf3_lista,
                           pathway.id = "hsa04966",
                           species    = "hsa",
                           out.suffix = "tm7sf3_Collecting_duct_acid_secretion")



