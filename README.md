# scSAGRN
Single-Cell Multi-Omics Gene Regulatory Network Analysis Pipeline
Overview
This pipeline is used to infer Gene Regulatory Networks (GRNs) from single-cell RNA-seq and ATAC-seq data, including peak-gene links and transcription factor (TF)-gene links. The pipeline is divided into the following steps:
1.Raw Data Preparation
  Paired gene expression data and chromatin accessibility data are used as input.
2.Data Preprocessing
  Refer to datapreprocess to preprocess raw single-cell multi-omics data and obtain highly variable genes.
3.Peak-Gene Link Inference
  Input processed seurat type data and highly variable genes, use runGenePeakLink function to obtain significant peak-gene relationships
  
  seurat <- readRDS(“A549_seurat_object.rds”)  
  genelist <- read.table(“highly_variable_genes.txt”)
  Peak_Gene <- runGenePeakLink(seurat.
                             genome = “hg38”, 
                             geneList = genelist[[1]],
                             windowPadSize = 150000,
                             nCores = 16,
                             p.cut = NULL, 
                             n_bg = 100)
4.Gene Score Calculation
  Perform Peak_Gene filtering based on p-value and run runGeneScore to get gene score
  
  Peak_Gene.filt <- Peak_Gene%>% filter(pvalZ <= 0.05)
  Gene <- unique(Peak_Gene.filt$Gene)
  GeneScore <- runGeneScore(seurat.                          
                          peakgene = Peak_Gene.filt,                         
                          geneList = Gene,                         
                          nCores = 8)

5.TF-Gene Link Inference
Input the seurat data and the peak-gene links obtained by filtering, and use the TFGeneGRN function to obtain the TF-gene relationship

TF_Gene <- TFGeneGRN(seurat, 
                     peakgene = Peak_Gene.filt, 
                     genome = “hg38”,
                     GeneScore = GeneScore,
                     nCores = 8)
