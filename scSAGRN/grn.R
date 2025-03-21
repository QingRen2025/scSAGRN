library(dplyr)
library(FNN)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)

RNA<-readRDS("method/pbmc/my/processdata/weighted_rna.rds")  
ATAC<- readRDS("method/pbmc/my/processdata/weighted_atac.SE.rds")  
dim(ATAC) # Peaks x ATAC cells
dim(RNA) # Genes x RNA cells
genelist <- read.table("method/pbmc/my/processdata/high_variable_genes.txt")

cisCorr <- runGenePeakLink(WNNATAC = ATAC,
                                 WNNRNA = RNA,
                                 genome = "hg38", # One of hg19, mm10 or hg38 
                                 geneList =genelist[[1]],
                                 windowPadSize = 150000,
                                 nCores = 8,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = 100)
head(cisCorr)
write.csv(cisCorr,file="method/pbmc/result/pbmc_peak_gene.csv",row.names=F)

cisCorr <- read.csv("method/pbmc/result/pbmc_peak_gene.csv")


cisCorr.filt <- cisCorr%>% filter(pvalZ <= 0.05)

dorcGenes <- GeneRankPlot(dorcTab = cisCorr.filt,
                        cutoff = 1, # No. sig peaks needed to be called a DORC
                        labelTop = 20,
                        returnGeneList = TRUE, # Set this to FALSE for just the plot
                        force=2)


# Unfiltered 获得基因名称
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

dorcMat <- getDORCScores(WNNATAC = ATAC,                          
                         dorcTab = cisCorr.filt,                         
                         geneList = dorcGenes,                         
                         nCores = 4)
dim(dorcMat)

saveRDS(dorcMat, file = "method/pbmc/result/dorcMat.rds")
# 使用readRDS函数加载RDS文件  
dorcMat <- readRDS(file = "method/pbmc/result/dorcMat.rds")

TF_Gene <- TFGeneGRN(WNNATAC = ATAC, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat,
                     WNNRNA = RNA, 
                     nCores = 4)
# 保存 TF_Gene 对象
saveRDS(TF_Gene, file = "method/pbmc/result/TF_Generesults.rds")
