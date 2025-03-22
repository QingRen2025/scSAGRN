runGenePeakLink <- function (seurat, genome, geneList, windowPadSize = 150000, 
                             normalizeATACmat = TRUE, nCores = 4, keepPosCorOnly = TRUE, 
                             keepMultiMappingPeaks = FALSE, n_bg = 100, p.cut = NULL) 
{
  
  if (!"RNA" %in% names(seurat@assays)) 
    stop("Seurat object must contain 'RNA' assay")
  if (!"ATAC" %in% names(seurat@assays)) 
    stop("Seurat object must contain 'ATAC' assay")
  if (!"wsnn" %in% names(seurat@graphs)) 
    stop("Seurat object must contain 'wsnn' graph in @graphs")
  
  wnn_matrix <- seurat@graphs$wsnn
  wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix)
  
  rna <- A549@assays[["RNA"]]@layers$data
  gene_names <- rownames(A549@assays[["RNA"]]@features@.Data)  
  rownames(rna) <- gene_names  #对应基准测试数据集的处理
  
  # RNA 加权平均
  weighted_rna <- rna %*% wnn_matrix_normalized
  # 查看结果的维度
  dim(weighted_rna)
  WNNRNA <- weighted_rna
  
  atac <- seurat@assays[["ATAC"]]@data
  # ATAC 加权平均
  weighted_atac <- atac %*% wnn_matrix_normalized
  # 查看结果的维度
  dim(weighted_atac)
  
  if (dim(weighted_rna)[2] != dim(weighted_atac)[2]) {
    stop("Error: The number of columns in 'weighted_rna' and 'weighted_atac' do not match.\n",
         "weighted_rna columns: ", dim(weighted_rna)[2], 
         "weighted_atac columns: ", dim(weighted_atac)[2])
  }
  
  peakRanges <- seurat@assays[["ATAC"]]@ranges
  peakCounts <- weighted_atac
  cellMeta <- DataFrame(seurat@meta.data)
  WNNATAC <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=peakCounts),
                                                        rowRanges = peakRanges,
                                                        colData = cellMeta)
  # 确保必要包已加载
  if (!requireNamespace("FigR", quietly = TRUE)) 
    stop("Please install FigR: devtools::install_github('author/FigR')")
  
  
  # 原始峰范围
  peakRanges.OG <- granges(WNNATAC)
  rownames(WNNATAC) <- paste0("Peak", 1:nrow(WNNATAC))
  ATACmat <- assay(WNNATAC)
  
  # 归一化 ATAC 数据
  if (normalizeATACmat) 
    ATACmat <- centerCounts(ATACmat)
  
  # RNA 基因名检查
  if (is.null(rownames(WNNRNA))) 
    stop("RNA matrix must have gene names as rownames")
  
  # 去除全零的峰和基因
  if (any(Matrix::rowSums(assay(WNNATAC)) == 0)) {
    message("Peaks with 0 accessibility across cells exist .. Removing them.")
    peaksToKeep <- Matrix::rowSums(assay(WNNATAC)) != 0
    WNNATAC <- WNNATAC[peaksToKeep, ]
    ATACmat <- ATACmat[peaksToKeep, ]
  }
  if (any(Matrix::rowSums(WNNRNA) == 0)) {
    message("Genes with 0 expression across cells exist .. Removing them.")
    genesToKeep <- Matrix::rowSums(WNNRNA) != 0
    WNNRNA <- WNNRNA[genesToKeep, ]
  }
  
  # 获取 TSS 注释
  if (!genome %in% c("hg38", "mm10")) 
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  
  switch(genome,  hg38 = {
    TSSg <- FigR::hg38TSSRanges
  }, mm10 = {
    TSSg <- FigR::mm10TSSRanges
  })
  
  names(TSSg) <- as.character(TSSg$gene_name)
  
  # 筛选 geneList 中有效基因
  if (!is.null(geneList)) {
    valid_genes <- geneList[geneList %in% names(TSSg)]
    if (length(valid_genes) == 0) 
      stop("None of the genes in the gene list are present in the TSS annotations.\n")
    TSSg <- TSSg[valid_genes]
  }
  
  # 基因交集
  genesToKeep <- intersect(names(TSSg), rownames(WNNRNA))
  WNNRNA <- WNNRNA[genesToKeep, ]
  TSSg <- TSSg[genesToKeep]
  
  # 保存有效基因列表
  youxiaogene <- genesToKeep
  
  # 构建 TSS 窗口
  TSSflank <- GenomicRanges::flank(TSSg, width = windowPadSize, both = TRUE)
  
  # 提取峰中心
  peakSummits <- resize(granges(WNNATAC), width = 1, fix = "center")
  
  # 基因与峰重叠
  genePeakOv <- findOverlaps(query = TSSflank, subject = peakSummits)
  
  # 保存有效峰列表
  youxiaopeak <- unique(subjectHits(genePeakOv))
  youxiaopeak_ranges <- peakRanges.OG[youxiaopeak]
  numPairs <- length(youxiaopeak)
  
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ", 
      length(genesToKeep), "\n")
  cat("Number of peaks overlapping gene TSS windows: ", 
      length(youxiaopeak), "\n")
  
  # 背景计算
  set.seed(123)
  cat("Determining background peaks ..\n")
  if (is.null(rowData(WNNATAC)$bias)) {
    if (genome %in% "mm10") 
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    WNNATAC <- chromVAR::addGCBias(WNNATAC, genome = myGenome)
  }
  
  cat("Using ", n_bg, " iterations ..\n\n")
  set.seed(123)
  bg <- chromVAR::getBackgroundPeaks(WNNATAC, niterations = n_bg)
  cat("Computing gene-peak correlations ..\n")
  pairsPerChunk <- 500
  largeChunkSize <- 5000
  startingPoint <- 1
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize - 1
  chunkEnds[length(chunkEnds)] <- numPairs
  library(doParallel)
  dorcList <- list()
  for (i in 1:length(chunkStarts)) {
    cat("Running pairs: ", chunkStarts[i], "to", chunkEnds[i], 
        "\n")
    ObsCor <- FigR::PeakGeneCor(ATAC = ATACmat, RNA = WNNRNA, 
                                OV = genePeakOv[chunkStarts[i]:chunkEnds[i]], chunkSize = pairsPerChunk, 
                                ncores = nCores, bg = bg)
    gc()
    dorcList[[i]] <- ObsCor
  }
  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)
  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg) + 3)
  if (keepPosCorOnly) {
    cat("Only considering positive correlations ..\n")
    dorcTab <- dorcTab %>% dplyr::filter(rObs > 0)
  }
  if (!keepMultiMappingPeaks) {
    cat("Keeping max correlation for multi-mapping peaks ..\n")
    dorcTab <- dorcTab %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs == 
                                                                     max(rObs))
  }
  dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]
  dorcTab$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTab$Peak], 
                                           "Peak", 2))
  dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[, 
                                                         permCols]))
  dorcTab$rBgMean <- rowMeans(dorcTab[, permCols])
  dorcTab$pvalZ <- 1 - stats::pnorm(q = dorcTab$rObs, mean = dorcTab$rBgMean, 
                                    sd = dorcTab$rBgSD)
  cat("\nFinished!\n")
  if (!is.null(p.cut)) {
    cat("Using significance cut-off of ", p.cut, " to subset to resulting associations\n")
    dorcTab <- dorcTab[dorcTab$pvalZ <= p.cut, ]
  }
  dorcTab$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[dorcTab$Peak])), 
                              paste(start(peakRanges.OG[dorcTab$Peak]), end(peakRanges.OG[dorcTab$Peak]), 
                                    sep = "-"), sep = ":")
  return(as.data.frame(dorcTab[, c("Peak", "PeakRanges", "Gene", 
                                   "rObs", "pvalZ")], stringsAsFactors = FALSE))
}
