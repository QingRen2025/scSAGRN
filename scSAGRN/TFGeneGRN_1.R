TFGeneGRN <-function (seurat, GeneK = 30, peakgene, n_bg = 50, genome, GeneScore, Genes = NULL, nCores = 1) 
{
  
  if (!"RNA" %in% names(seurat@assays)) 
    stop("Seurat object must contain 'RNA' assay")
  if (!"ATAC" %in% names(seurat@assays)) 
    stop("Seurat object must contain 'ATAC' assay")
  if (!"wsnn" %in% names(seurat@graphs)) 
    stop("Seurat object must contain 'wsnn' graph in @graphs")
  
  wnn_matrix <- seurat@graphs$wsnn
  wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix)
  
  rna <- seurat@assays[["RNA"]]@layers$data
  gene_names <- rownames(seurat@assays[["RNA"]]@features@.Data)  
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
  stopifnot(all.equal(ncol(GeneScore), ncol(WNNRNA)))
  if (!all(c("Peak", "Gene") %in% colnames(peakgene))) 
    stop("Expecting fields Peak and Gene in peakgene data.frame .. see runGenePeakcorr function in BuenRTools")
  if (all(grepl("chr", peakgene$Peak, ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")
    if (!(all(grepl("chr", rownames(WNNATAC), ignore.case = TRUE)))) 
      stop("Peak regions provided in peakgene data.frame but not found as rownames in input SE")
    if (!all(peakgene$Peak %in% rownames(WNNATAC))) 
      stop("Found Gene peak region not present in input SE.. make sure Gene calling output corresponds to same input SE as the one provided here ..")
  }
  else {
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
    if (max(peakgene$Peak) > nrow(WNNATAC)) 
      stop("Found Gene peak index outside range of input SE.. make sure Gene calling output corresponds to same input SE as the one provided here ..")
  }
  if (is.null(Genes)) {
    Genes <- rownames(GeneScore)
  }
  else {
    cat("Using specified list of genes ..\n")
    if (!(all(Genes %in% rownames(GeneScore)))) {
      cat("One or more of the gene names supplied is not present in the Gene matrix provided: \n")
      cat(Genes[!Genes %in% rownames(GeneScore)], 
          sep = ", ")
      cat("\n")
      stop()
    }
  }
  Gene.knn <- FNN::get.knn(data = t(scale(Matrix::t(GeneScore))), 
                           k = GeneK)$nn.index
  rownames(Gene.knn) <- rownames(GeneScore)
  if (is.null(SummarizedExperiment::rowData(WNNATAC)$bias)) {
    if (genome %in% "mm10") 
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    if (genome %in% "hg19") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    WNNATAC <- chromVAR::addGCBias(WNNATAC, genome = myGenome)
  }
  packagePath <- find.package("FigR", lib.loc = NULL, quiet = TRUE)
  if (grepl("hg", genome)) {
    pwm <- readRDS(paste0(packagePath, "/data/cisBP_human_pfms_2021.rds"))
  }
  else {
    pwm <- readRDS(paste0(packagePath, "/data/cisBP_mouse_pfms_2021.rds"))
  }
  if (all(grepl("_", names(pwm), fixed = TRUE))) 
    names(pwm) <- FigR::extractTFNames(names(pwm))
  message("Removing genes with 0 expression across cells ..\n")
  WNNRNA <- WNNRNA[Matrix::rowSums(WNNRNA) != 0, ]
  myGeneNames <- gsub(x = rownames(WNNRNA), pattern = "-", 
                      replacement = "")
  rownames(WNNRNA) <- myGeneNames
  motifsToKeep <- intersect(names(pwm), myGeneNames)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = WNNATAC, pwms = pwm[motifsToKeep], 
                                       genome = genome)
  motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix)) != 
                         0]
  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if (any(Matrix::rowSums(assay(WNNATAC)) == 0)) {
    ATAC.mat <- assay(WNNATAC)
    ATAC.mat <- cbind(ATAC.mat, 1)
    WNNATAC.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = ATAC.mat), 
                                                              rowRanges = granges(WNNATAC))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(WNNATAC.new, niterations = n_bg)
  }
  else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(WNNATAC, niterations = n_bg)
  }
  cat("Testing ", length(motifsToKeep), " TFs\n")
  cat("Testing ", nrow(GeneScore), " Genes\n")
  library(doParallel)
  if (nCores > 1) 
    message("Running FigR using ", nCores, " cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(Genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  mZtest.list <- foreach(g = Genes, .options.snow = opts, 
                         .packages = c("FigR", "dplyr", "Matrix", "Rmpfr")) %dopar% 
    {
      GeneNNpeaks <- unique(peakgene$Peak[peakgene$Gene %in% 
                                           c(g, rownames(GeneScore)[Gene.knn[g, ]])])
      if (usePeakNames) 
        GeneNNpeaks <- which(rownames(WNNATAC) %in% GeneNNpeaks)
      mZ <- FigR::motifPeakZtest(peakSet = GeneNNpeaks, 
                                 bgPeaks = bg, tfMat = assay(motif_ix))
      mZ <- mZ[, c("gene", "z_test")]
      colnames(mZ)[1] <- "Motif"
      colnames(mZ)[2] <- "Enrichment.Z"
      mZ$Enrichment.P <- 2 * pnorm(abs(mZ$Enrichment.Z), 
                                   lower.tail = FALSE)
      mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
      mZ <- cbind(Gene = g, mZ)
      corr.r <- cor(GeneScore[g, ], t(as.matrix(WNNRNA[mZ$Motif, 
      ])), method = "spearman")
      stopifnot(all.equal(colnames(corr.r), mZ$Motif))
      mZ$Corr <- corr.r[1, ]
      mZ$Corr.Z <- scale(mZ$Corr, center = TRUE, scale = TRUE)[, 
                                                               1]
      mZ$Corr.P <- 2 * pnorm(abs(mZ$Corr.Z), lower.tail = FALSE)
      mZ$Corr.log10P <- sign(mZ$Corr.Z) * -log10(mZ$Corr.P)
      return(mZ)
    }
  cat("Finished!\n")
  cat("Merging results ..\n")
  TFenrich.d <- do.call("rbind", mZtest.list)
  dim(TFenrich.d)
  rownames(TFenrich.d) <- NULL
  TFenrich.d <- TFenrich.d %>% dplyr::mutate(Score = sign(Corr) * 
                                               as.numeric(-log10(1 - (1 - Rmpfr::mpfr(Enrichment.P, 
                                                                                      100)) * (1 - Rmpfr::mpfr(Corr.P, 100)))))
  TFenrich.d$Score[TFenrich.d$Enrichment.Z < 0] <- 0
  colnames(TFenrich.d)[colnames(TFenrich.d) == "DORC"] <- "Gene"
  colnames(TFenrich.d)[colnames(TFenrich.d) == "Motif"] <- "TF"
  TFenrich.d
  
}
