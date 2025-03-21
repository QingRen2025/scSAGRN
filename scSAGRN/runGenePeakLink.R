runGenePeakLink <- function (WNNATAC, WNNRNA, genome, geneList = NULL, windowPadSize = 150000, 
                             normalizeATACmat = TRUE, nCores = 4, keepPosCorOnly = TRUE, 
                             keepMultiMappingPeaks = FALSE, n_bg = 100, p.cut = NULL) 
{
  stopifnot(inherits(WNNATAC, "RangedSummarizedExperiment"))
  stopifnot(inherits(WNNRNA, c("Matrix", "matrix")))
  if (!all.equal(ncol(WNNATAC), ncol(WNNRNA))) 
    stop("Input ATAC and RNA objects must have same number of cells")
  message("Assuming paired scATAC/scRNA-seq data ..")
  
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
  if (!genome %in% c("hg19", "hg38", "mm10")) 
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  
  switch(genome, hg19 = {
    TSSg <- FigR::hg19TSSRanges
  }, hg38 = {
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
  write.table(youxiaogene, file = "method/pbmc/my/processdata/youxiaogene.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 构建 TSS 窗口
  TSSflank <- GenomicRanges::flank(TSSg, width = windowPadSize, both = TRUE)
  
  # 提取峰中心
  peakSummits <- resize(granges(WNNATAC), width = 1, fix = "center")
  
  # 基因与峰重叠
  genePeakOv <- findOverlaps(query = TSSflank, subject = peakSummits)
  
  # 保存有效峰列表
  youxiaopeak <- unique(subjectHits(genePeakOv))
  youxiaopeak_ranges <- peakRanges.OG[youxiaopeak]
  write.table(as.data.frame(youxiaopeak_ranges), file = "method/pbmc/my/processdata/youxiaopeak.txt", row.names = FALSE, quote = FALSE, sep = "\t")
  
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ", 
      length(genesToKeep), "\n")
  cat("Number of peaks overlapping gene TSS windows: ", 
      length(youxiaopeak), "\n")
  
  # 背景计算
  set.seed(123)
  cat("Determining background peaks ..\n")
  if (is.null(rowData(WNNATAC)$bias)) {
    if (genome %in% "hg19") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
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
