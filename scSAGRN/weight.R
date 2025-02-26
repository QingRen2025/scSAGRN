# 加载必要的库
library(SummarizedExperiment)
library(Matrix)
library(GenomicRanges)

wnn_matrix <- readRDS("method/pbmc/my/processdata/wnn_matrix.rds")
head(wnn_matrix)
# 对每一行进行标准化（使每一行的和为 1）
wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix)
# 保存标准化后的 WNN 矩阵为 RDS 文件
saveRDS(wnn_matrix_normalized, file = "method/pbmc/my/processdata/wnn_matrix_normalized.rds")
# 加载 RDS 文件
wnn_matrix_normalized <- readRDS("method/pbmc/my/processdata/wnn_matrix_normalized.rds")

pbmc <- readRDS("method/pbmc/my/processdata/pbmc_seurat_object.rds")

rna <- pbmc@assays$SCT@data
# RNA 加权平均
weighted_rna <- rna %*% wnn_matrix_normalized
# 查看结果的维度
dim(weighted_rna)
# 保存 weighted_rna 为 RDS 文件
saveRDS(weighted_rna, file = "method/pbmc/my/processdata/weighted_rna.rds")

atac <- pbmc@assays[["ATAC"]]@data
# ATAC 加权平均
weighted_atac <- atac %*% wnn_matrix_normalized
# 查看结果的维度
dim(weighted_atac)
saveRDS(weighted_atac, file = "method/pbmc/my/processdata/weighted_atac.rds")

peakRanges <- pbmc@assays[["ATAC"]]@ranges

peakCounts <- weighted_atac


cellMeta <- DataFrame(pbmc@meta.data)


weighted_atac.SE <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=peakCounts),
                                                 rowRanges = peakRanges,
                                                 colData = cellMeta)
saveRDS(weighted_atac.SE, file = "method/pbmc/my/processdata/weighted_atac.SE.rds")


