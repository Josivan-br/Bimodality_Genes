library(stringr)
library(stringr)
library(data.table)
library(dplyr)

##### Sobreposição P2

if(!require('UpSetR')) { install.packages('UpSetR') }

setwd("E:/AmostrasKIRC/KIRC_Script/AlgotitmoAmostras/samples/KIRC/P2")

files <- list.files(pattern = "*_0.02_2.lst")

genes <- str_replace(files, "_0.02_2.lst", "")

lista <- data.frame() # cria um data.frame vazio para todos as listas de cada gene

for(i in files) {
  temp <- fread(i, header = F)# ler cada arq. em um arq. temporario
  #separa o nome do arquivo strsplit("-") e guarda para cada arquivo i 
  temp$gene <- strsplit(i, split = "*_0.02_2.lst")[[1]][1] 
  
  lista <- rbind(lista, temp) # sobrepoe os arq. temp data.frame lista
}
### Mudar o numero de colunas
up_set <- data.frame(matrix(nrow=length(unique(lista$V1)), ncol= 6, dimnames = list(unique(lista$V1), genes) )) # cria um data.frame vazio
up_set[, genes] <- 0

for(g in genes){
 list_g  <- read.table(paste0(g, "_0.02_2.lst"), header = F, stringsAsFactors = F)
 up_set[list_g$V1, g] <- 1
}

upset(up_set, sets = genes, sets.bar.color = "#56B4E9",  order.by = c("freq", "degree"),
      text.scale = 1.2,  set_size.show = TRUE, set_size.scale_max = 450,
      mainbar.y.label = "Intersection Size of Unique Patiens")

png(file="UpSet_5Genes.png", width=800, height=800)
print(
  upset(up_set, sets = genes, sets.bar.color = "#56B4E9",  order.by = c("freq", "degree"),
        text.scale = 2,  set_size.show = TRUE, set_size.scale_max = 450,
        mainbar.y.label = "Intersection Size of Unique Patiens")
      )

dev.off()


##### Sobreposição P1

if(!require('UpSetR')) { install.packages('UpSetR') }

setwd("E:/AmostrasKIRC/KIRC_Script/AlgotitmoAmostras/samples")

files <- list.files(pattern = "*_0.02_1.lst")

genes <- str_replace(files, "_0.02_1.lst", "")

lista <- data.frame() # cria um data.frame vazio para todos as listas de cada gene

for(i in files) {
  temp <- fread(i, header = F)# ler cada arq. em um arq. temporario
  #separa o nome do arquivo strsplit("-") e guarda para cada arquivo i 
  temp$gene <- strsplit(i, split = "*_0.02_1.lst")[[1]][1] 
  
  lista <- rbind(lista, temp) # sobrepoe os arq. temp data.frame lista
}
### Mudar o numero de colunas
up_set <- data.frame(matrix(nrow=length(unique(lista$V1)), ncol= 6, dimnames = list(unique(lista$V1), genes) )) # cria um data.frame vazio
up_set[, genes] <- 0

for(g in genes){
  list_g  <- read.table(paste0(g, "_0.02_1.lst"), header = F, stringsAsFactors = F)
  up_set[list_g$V1, g] <- 1
}

upset(up_set, sets = genes, sets.bar.color = "#56B4E9",  order.by = c("freq", "degree"),
      text.scale = 1.2,  set_size.show = TRUE, set_size.scale_max = 450,
      mainbar.y.label = "Intersection Size of Unique Patiens")

png(file="UpSet_5Genes.png", width=800, height=800)
print(
  upset(up_set, sets = genes, sets.bar.color = "#56B4E9",  order.by = c("freq", "degree"),
        text.scale = 2,  set_size.show = TRUE, set_size.scale_max = 450,
        mainbar.y.label = "Intersection Size of Unique Patiens")
)

dev.off()


