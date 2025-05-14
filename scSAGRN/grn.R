library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)

pbmc <- readRDS("method/pbmc/processdata/pbmc__seurat_object.rds")  
genelist <- read.table("method/pbmc/processdata/high_variable_genes.txt")

Peak_Gene <- runGenePeakLink(pbmc,
                           genome = "hg38", 
                           geneList =genelist[[1]],
                           windowPadSize = 150000,
                           nCores = 8,
                           p.cut = NULL, 
                           n_bg = 100)
head(Peak_Gene)
write.csv(Peak_Gene,file="method/pbmc/result/pbmc_peak_gene.csv",row.names=F)

#根据P值进行过滤
Peak_Gene.filt <- Peak_Gene%>% filter(pvalZ <= 0.05)
Gene <- unique(Peak_Gene.filt$Gene)

GeneScore <- runGeneScore(pbmc,                          
                         peakgene = Peak_Gene.filt,                         
                         geneList = Gene,                         
                         nCores = 4)

saveRDS(GeneScore, file = "method/pbmc/result/GeneScore.rds")

TF_Gene <- TFGeneGRN(pbmc, 
                     peakgene = Peak_Gene.filt, 
                     genome = "hg38",
                     GeneScore = GeneScore,
                     nCores = 4)

saveRDS(TF_Gene, file = "method/pbmc/result/TF_Generesults.rds")
