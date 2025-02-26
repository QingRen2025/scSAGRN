TFGeneGRN <-function (ATAC.se, dorcK = 30, dorcTab, n_bg = 50, genome, dorcMat, 
          rnaMat, dorcGenes = NULL, nCores = 1) 
{
  stopifnot(all.equal(ncol(dorcMat), ncol(rnaMat)))
  if (!all(c("Peak", "Gene") %in% colnames(dorcTab))) 
    stop("Expecting fields Peak and Gene in dorcTab data.frame .. see runGenePeakcorr function in BuenRTools")
  if (all(grepl("chr", dorcTab$Peak, ignore.case = TRUE))) {
    usePeakNames <- TRUE
    message("Detected peak region names in Peak field")
    if (!(all(grepl("chr", rownames(ATAC.se), ignore.case = TRUE)))) 
      stop("Peak regions provided in dorcTab data.frame but not found as rownames in input SE")
    if (!all(dorcTab$Peak %in% rownames(ATAC.se))) 
      stop("Found DORC peak region not present in input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }
  else {
    usePeakNames <- FALSE
    message("Assuming peak indices in Peak field")
    if (max(dorcTab$Peak) > nrow(ATAC.se)) 
      stop("Found DORC peak index outside range of input SE.. make sure DORC calling output corresponds to same input SE as the one provided here ..")
  }
  if (is.null(dorcGenes)) {
    dorcGenes <- rownames(dorcMat)
  }
  else {
    cat("Using specified list of dorc genes ..\n")
    if (!(all(dorcGenes %in% rownames(dorcMat)))) {
      cat("One or more of the gene names supplied is not present in the DORC matrix provided: \n")
      cat(dorcGenes[!dorcGenes %in% rownames(dorcMat)], 
          sep = ", ")
      cat("\n")
      stop()
    }
  }
  DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))), 
                           k = dorcK)$nn.index
  rownames(DORC.knn) <- rownames(dorcMat)
  if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
    if (genome %in% "hg19") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if (genome %in% "mm10") 
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38") 
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
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
  rnaMat <- rnaMat[Matrix::rowSums(rnaMat) != 0, ]
  myGeneNames <- gsub(x = rownames(rnaMat), pattern = "-", 
                      replacement = "")
  rownames(rnaMat) <- myGeneNames
  motifsToKeep <- intersect(names(pwm), myGeneNames)
  cat("Getting peak x motif matches ..\n")
  motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se, pwms = pwm[motifsToKeep], 
                                       genome = genome)
  motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix)) != 
                         0]
  cat("Determining background peaks ..\n")
  cat("Using ", n_bg, " iterations ..\n\n")
  if (any(Matrix::rowSums(assay(ATAC.se)) == 0)) {
    ATAC.mat <- assay(ATAC.se)
    ATAC.mat <- cbind(ATAC.mat, 1)
    ATAC.se.new <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = ATAC.mat), 
                                                              rowRanges = granges(ATAC.se))
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se.new, niterations = n_bg)
  }
  else {
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
  }
  cat("Testing ", length(motifsToKeep), " TFs\n")
  cat("Testing ", nrow(dorcMat), " DORCs\n")
  library(doParallel)
  if (nCores > 1) 
    message("Running FigR using ", nCores, " cores ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(nCores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  mZtest.list <- foreach(g = dorcGenes, .options.snow = opts, 
                         .packages = c("FigR", "dplyr", "Matrix", "Rmpfr")) %dopar% 
    {
      DORCNNpeaks <- unique(dorcTab$Peak[dorcTab$Gene %in% 
                                           c(g, rownames(dorcMat)[DORC.knn[g, ]])])
      if (usePeakNames) 
        DORCNNpeaks <- which(rownames(ATAC.se) %in% DORCNNpeaks)
      mZ <- FigR::motifPeakZtest(peakSet = DORCNNpeaks, 
                                 bgPeaks = bg, tfMat = assay(motif_ix))
      mZ <- mZ[, c("gene", "z_test")]
      colnames(mZ)[1] <- "Motif"
      colnames(mZ)[2] <- "Enrichment.Z"
      mZ$Enrichment.P <- 2 * pnorm(abs(mZ$Enrichment.Z), 
                                   lower.tail = FALSE)
      mZ$Enrichment.log10P <- sign(mZ$Enrichment.Z) * -log10(mZ$Enrichment.P)
      mZ <- cbind(DORC = g, mZ)
      corr.r <- cor(dorcMat[g, ], t(as.matrix(rnaMat[mZ$Motif, 
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
