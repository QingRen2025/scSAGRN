runGeneScore <- function (seurat, peakgene, normalizeATACmat = TRUE, geneList, 
                          nCores = 4) 
{
  if (!"ATAC" %in% names(seurat@assays)) 
    stop("Seurat object must contain 'ATAC' assay")
  if (!"wsnn" %in% names(seurat@graphs)) 
    stop("Seurat object must contain 'wsnn' graph in @graphs")
  
  wnn_matrix <- seurat@graphs$wsnn
  wnn_matrix_normalized <- wnn_matrix / rowSums(wnn_matrix)
  
  atac <- seurat@assays[["ATAC"]]@data
  # ATAC 加权平均
  weighted_atac <- atac %*% wnn_matrix_normalized
  
  peakRanges <- seurat@assays[["ATAC"]]@ranges
  peakCounts <- weighted_atac
  cellMeta <- DataFrame(seurat@meta.data)
  WNNATAC <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=peakCounts),
                                                        rowRanges = peakRanges,
                                                        colData = cellMeta)
  # 确保必要包已加载
  if (!requireNamespace("FigR", quietly = TRUE)) 
    stop("Please install FigR: devtools::install_github('author/FigR')")
  if (!all(c("Peak", "Gene") %in% colnames(peakgene))) 
    stop("The provided gene-peak table must have columns named Peak and Gene ..")
  if (any(peakgene$Peak > nrow(WNNATAC))) 
    stop("One or more peak indices in the gene-peak table are larger than the total number of peaks in the provided ATAC SE object ..\n Make sure the exact same SummarizedExperiment object is provided here as was for running the runGenePeakcorr function ..\n")
  if (!is.null(geneList)) {
    if (!(all(geneList %in% as.character(peakgene$Gene)))) 
      stop("One or more of the gene names supplied is not present in the gene-peak table provided..\n")
    if (length(geneList) > 50) {
      message("Running gene scoring for ", length(geneList), 
              " genes: ", paste(geneList[1:20], collapse = ", "), 
              ", ... , ... , ... (truncated display)")
    }
    else {
      message("Running gene scoring for ", length(geneList), 
              " genes: ", paste(geneList, collapse = "\n"))
    }
    cat("........\n")
    peakgene <- peakgene[peakgene$Gene %in% geneList, ]
    Genes <- sort(as.character(unique(peakgene$Gene)))
  }
  else {
    Genes <- sort(as.character(unique(peakgene$Gene)))
    cat("Running gene scoring for all genes in annotation! (n = ", 
        length(Genes), ")\n", sep = "")
  }
  if (normalizeATACmat) {
    cat("Normalizing scATAC counts ..\n")
    ATAC.mat <- assay(centerCounts(WNNATAC, chunkSize = 5000))
    gc()
  }
  else {
    cat("Assuming provided scATAC counts are normalized ..\n")
    ATAC.mat <- assay(WNNATAC)
  }
  time_elapsed <- Sys.time()
  if (Sys.info()["sysname"] %in% "Windows") {
    message("Windows OS detected .. Cannot support parallilzation using mclapply for mc.cores > 1")
    message("Using 1 core instead ..\n")
    nCores <- 1
  }
  cat("Computing gene scores ..\n")
  cat("Running in parallel using ", nCores, "cores ..\n")
  geneMatL <- pbmcapply::pbmclapply(X = Genes, FUN = function(x) {
    genePeaks <- unique(peakgene$Peak[peakgene$Gene %in% x])
    if (length(genePeaks) > 1) {
      geneCounts <- Matrix::colSums(ATAC.mat[genePeaks, 
      ])
    }
    else if (length(genePeaks == 1)) {
      geneCounts <- ATAC.mat[genePeaks, ]
    }
  }, mc.cores = nCores)
  geneMat <- Matrix::Matrix(do.call("rbind", geneMatL), sparse = TRUE)
  rownames(geneMat) <- Genes
  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)), 
      "\n\n")
  return(geneMat)
}