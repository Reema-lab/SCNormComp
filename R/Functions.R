#' A test package
#' 
#' 
# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Fahrenheit conversion
#'
#' Convert degrees Fahrenheit temperatures to degrees Celsius
#' @param Exp_mar The temperature in degrees Fahrenheit
#' @return The expression matrix after imputation
#' @examples 
#' temp1 <- normRes(50);
#' temp2 <- compRes(70) );
#' temp3 <- impRes(65)
#' 
#' @export
impRes <- function(x,y) {
  #Create a matrix of C rows and M columns, to store the genes corresponding to each cluster.
  #Get the Assay Data from the Seurat object.
  expression_data <- x
  n_row <- y
  #n_row <- length(unique(seuratobj_data$seurat_clusters))
  n_row <- y
  n_col <- ncol(expression_data)
  n_r <- nrow(expression_data)
  
  matlist <- matrix(nrow=n_row,ncol = n_col)
  ind <- 1
  nc <- n_row-1
  for(i in 0:nc)
  {
    x <- which(cluster_nos==i)
    print(i)
    lx <- length(x)
    matlist[ind,1:lx] <- x 
    ind <- ind+1
  }
  
  expmat_copy2 <- matrix(nrow = n_r,ncol = n_col)
  expmat_copy <- matrix(nrow = n_r,ncol = n_col)
  expmat_copy2 <- expression_data
  
  for(g in 1:n_r)
  {
    prevlen <- 0
    
    for(i in 1:cluster_nos)
    {
      genelist <- c()
      a <- na.omit(matlist[i,])
      alen <- length(a)
      newlen <- prevlen + alen
      for(j in 1:alen)
        genelist <- append(genelist,expmat_copy2[g,a[j]])   
      if(i==1)
        expmat_copy2[g,1:newlen] <- genelist
      else  
        expmat_copy2[g,1:newlen] <- c(expmat_copy2[g,1:prevlen],genelist)
      prevlen <- newlen
    }
    print(g)
  }  
  
  #write.csv(expmat_copy2,"Expression Matrix Copy 2.")
  #exp_cp2 <- read.csv("Expression Matrix Copy 2.")
  #Run the different normalisation methods.
  
  copyexp2 <- expmat_copy2
  #Method: Proposed, ProbSCNorm
  #Extract the set of zero counts, N_{cz}
  
  #1. Store the number of cells of each cluster in a vector
  clus_len <- c()
  for(i in 1:cluster_nos)
  {
    clus_len[i] <- length(na.omit(matlist[i,]))
    
  }
  
  #Plot mean vs. standard deviation, before and after.
  #library(vsn)
  #meanSdPlot(as.matrix(expmat_copy2,T,"Mean","Standard Deviation")
  
  
  
  #nr <- nrow(expmat_copy2)
  #nc <- ncol(expmat_copy2)
  rowdropped <- c()
  #Store a copy of expmat_copy2 in copyexp2
  copyexp2 <- expmat_copy2
  #Retain genes with non-zero values in $75%$ of the cells
  for(i in 1:nr)
  {
    #count the number of zeroes
    cz <- sum(expmat_copy2[i,]==0)
    ncz <- sum(expmat_copy2[i,]!=0)
    
    ratio <- ncz/nc
    if(ratio<0.02)
    {  
      expmat_copy2 <- expmat_copy2[-i,]
      rowdropped <- append(i,rowdropped)
    }  
  }
  
  #meanSdPlot(as.matrix(copyexp2),T,"Mean","Standard Deviation")
  library(gamlss)
  F_fam <- c()
  #expmat_copy2 <- expression_data_unnormalized[1:1000,]
  #Now estimate the counts using the Binomial and Poisson distributions
  for(i in 1:n_r)
  {
    zeroind_list <- which(expmat_copy2[i,]==0)
    #Extract values in the index set of zero counts
    zero_c <- expmat_copy2[i,zeroind_list]
    
    nonzeroind_list <- which(expmat_copy2[i,]!=0)
    nonzero_c <- expmat_copy2[i,nonzeroind_list]
    l <- length(zeroind_list)
    
    
    #Run a binomial distribution on the set of zero counts.
    binom_prob <- dom_extr(expression_data)
    binom_res <- dbinom(zero_c,length(zero_c),binom_prob)
    
    prob_nonzero <- which(binom_res==1)
    
    #For those that have a non-zero count distribution, fit them into the other set using a Normal Distribution
    
    expmat_copy2
    
    F <- fitDist(nonzero_c, k = 2, type = "counts")
    F_fam[i] <- F$family[[1]]
    
    lzc <- length(zero_c)
    lnzc <- length(nonzero_c)
    lpnz <- length(prob_nonzero)
    print(paste0("Began",i))
    
    result = switch(   
      F_fam[i],   
      "LG" = { 
        miss_val <- rlnorm(lnzc)
      },   
      "ZIPF"= {
        miss_val <- rZIPF(lnzc)},   
      "GEOM"= {
        set.seed(10)
        miss_val <- rgeom(lnzc,0.5)},
      "ZAP" = {
        miss_val <- rZAP(lnzc, mu = 5, sigma = 0.1)
      },
      "ZIP2" = {
        miss_val <- rZIP2(lnzc, mu = 5, sigma = 0.1)
      },
      "PO" =   {
        miss_val <- rpo(lnzc,1)},
      "NBI" = {
        miss_val <- rNBI(lnzc,1,1)},
      "GPO" = {
        miss_val <- rGPO(lnzc,1,1)},
      "PIG" = {
        miss_val <- rPIG(lnzc,1,1)},
      "ZALG" = {
        miss_val <- rZAIG(lnzc, mu = 1, sigma = 1)},
      "ZINBI" = {
        miss_val <- rZINBI(lnzc, mu = 1, sigma = 1, nu = 0.5)},
      "ZANBI" = {
        miss_val <- rZANBI(lnzc, mu = 1, sigma = 1, nu = 0.5)},
      "ZIPIG" = {
        miss_val <- rZIPIG(lnzc, mu = 1, sigma = 1, nu = 0.5)},
      "ZAPIG " = {
        print("Here")
        miss_val <- rZAPIG(lnzc, mu = 1, sigma = 1, nu = 0.5)
      },
      "NBF" = {
        miss_val <- rNBF(lnzc, mu = 1, sigma = 1, nu = 0.5)
      },
      " SI" = {
        miss_val <- rSI(lnzc, mu = 0.5, sigma = 0.02, nu = -0.5)},
      "BNB" = {
        miss_val <- rBNB(lnzc, mu = 1, sigma = 1, nu = 1)
      },
      "ZINBF" = {
        miss_val <- rZINBF(lnzc, mu = 1, sigma = 1, nu = 2, tau = 0.1)    
      },
      "ZASICHEL" = {
        miss_val <- rZASICHEL(lnzc, mu = 1, sigma = 1, nu = -0.5, tau = 0.1)
      },
      "ZISICHEL " = {
        miss_val <- rSICHEL(lnzc, mu=1, sigma=1, nu=-0.5)
      },
      "ZABNB " = {
        miss_val <- rZABNB(lnzc, mu = 1, sigma = 1, nu = 1, tau = 0.1)
      },
      "ZIBNB" = {
        miss_val <- rZIBNB(lnzc, mu = 1, sigma = 1, nu = 1, tau = 0.1)  
      },
      "ZAZIPF" = {
        miss_val <-  rZAZIPF(lnzc, mu = 0.5, sigma = 0.1)
      }
    )   
    nonzero_c <- miss_val
    expmat_copy2[i,nonzeroind_list] <- miss_val
    print(paste0("Ended",i))
    
  }
  return(expmat_copy2)
}

#' Celsius conversion
#'
#' Compress a given gene expression matrix
#' @param X The expression matrix
#' @return The compressed matrix
#' @examples 
#' temp1 <- compRes(22);

#' @export
compRes <- function(x)
{
  cc <- c(length=nrow(x))
  nr <- nrow(x)
  nc <- ncol(x)
  for(i in 1:nrow(x))
  {
    for(j in  1:ncol(x))
    {
      #Calc the compression threshold
      zero_ind <- which(x[i,j]==0)
      zl <- length(zero_ind)
      cc <- zl/nc
      thr <- cc*median(x[i,])
      str_zero <- paste(zero_ind,collapse = ";")
      
    }
  }
}


#' Celsius conversion
#'
#' Find the cell types from any expression matrix
#' @param X The expression matrix
#' @return The compressed matrix
#' @examples 
#' temp1 <- compRes(22);

#' @export

domcell_types <- function(x)
{
  # script to annotate cell types from 20k Human PBMCs from a healthy female donor
  # setwd("~/Desktop/demo/singleCell_singleR_part2/scripts")
  
  library(SingleR)
  library(celldex)
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
  
  # 10X CellRanger .HDF5 format ---------
  hdf5_obj <- Read10X_h5(filename = "C:\\Users\\user\\Downloads\\20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5",
                         use.names = TRUE,
                         unique.features = TRUE)
  pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
  pbmc.seurat
  
  # QC and Filtering ----------
  # explore QC
  pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
 # VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
  pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 & mitoPercent < 10)
  
  # It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.
  
  # pre-process standard workflow ---------------
  pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
  #ElbowPlot(pbmc.seurat.filtered)
  pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
  pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
  pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)
  
  # running steps above to get clusters
  #DimPlot(pbmc.seurat.filtered, reduction = "umap")
  #View(pbmc.seurat.filtered@meta.data)
  
  # run SingleR with multiple reference datasets (default mode) ---------
  
  # for pbmc data, we will use two datasets
  hpca <- celldex::HumanPrimaryCellAtlasData()
  dice <- celldex::DatabaseImmuneCellExpressionData()
  
  # ...1. Strategy 1: Using reference-specific labels ----------
  hpca$label.main
  dice$label.main
  
  # adding ref info to labels
  hpca$label.main <- paste0('HPCA.', hpca$label.main)
  dice$label.main <- paste0('DICE.', dice$label.main)
  
  # create a combined ref based on shared genes
  shared <- intersect(rownames(hpca), rownames(dice))
  combined <- cbind(hpca[shared,], dice[shared,])
  combined
  combined$label.main
  
  # run singleR using combined ref
  # savings counts into a separate object
  pbmc_counts <- GetAssayData(pbmc.seurat.filtered, layer = 'counts')
  
  com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
  table(com.res1$labels)
  
  pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res1)), 'labels']
  #View(pbmc.seurat.filtered@meta.data)
  
  #DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)
  
  # ...2. Strategy 2: Comparing scores across references ----------
  
  hpca$label.main
  dice$label.main
  hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
  dice$label.main <- gsub('DICE\\.','', dice$label.main)
  
  com.res2 <- SingleR(test = pbmc_counts, 
                      ref = list(HPCA = hpca, DICE = dice),
                      labels = list(hpca$label.main, dice$label.main))
  
  # Check the final label from the combined assignment.
  table(com.res2$labels)
  
  # which reference scored best for which label?
  grouping <- paste0(com.res2$labels,'.', com.res2$reference)
  best_ref <- as.data.frame(split(com.res2, grouping))
  
  # get de. genes from each individual references
  metadata(com.res2$orig.results$HPCA)$de.genes
  metadata(com.res2$orig.results$DICE)$de.genes
  
  library(viridis)
  # Combined diagnostics
  #plotScoreHeatmap(com.res2)
  
  
  # ...3. Strategy 3: Using Harmonized Labels ----------
  
  hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
  dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')
  
  # Using the same sets of genes:
  shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
  hpca.ont <- hpca.ont[shared,]
  dice.ont <- dice.ont[shared,]
  
  # Showing the top 10 most frequent terms:
  tail(sort(table(hpca.ont$label.ont)),10)
  tail(sort(table(dice.ont$label.ont)), 10)
  
  # using label.ont instead on label.main while running SingleR
  
  com.res3 <- SingleR(test = pbmc_counts,
                      ref = list(HPCA = hpca.ont, DICE = dice.ont),
                      labels = list(hpca.ont$label.ont, dice.ont$label.ont))
  
  table(com.res3$labels)
  
  # How to map cell ontology terms? ----------------
  
  colData(hpca.ont)
  colData(dice.ont)
  
  hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
  hpca.mapping <- read.delim(hpca.fle, header = F)
  return(com.res3$labels)
}


#' Celsius conversion
#'
#' Find the cell types from any expression matrix
#' @param X The expression matrix
#' @return Sigmoid matrix
#'  @examples 
#' temp1 <- compRes(22);

#' @export

domextr <- function(x,y)
{
  expression_data <- x
  species <- y
  #Install the ExpressionAtlas package
  BiocManager::install("ExpressionAtlas")
  library(ExpressionAtlas)
  supressMessages(library(ExpressionAtlas))
  
  
  #Perform normalization using domain knowledge.
  BiocManager::install("ExpressionAtlas")
  library(glmGamPoi)
  library(ExpressionAtlas)
  
  #Set an evaluated parameter based upon the type of cell
  #Create a parameter vector to hold the parameter values for all the 19 possibilities (special cases: 4,10,18,19).
  #param <- c("Mesophyll","Mesophyll","Mesophyll","BS,XP1,XP2","Mesophyll","Mesophyll","Mesophyll","Mesophyll","Mesophyll","PP1,PC-XP,PC-PP","Mesophyll","Mesophyll","Epidermis","Mesophyll","Companion","Guard cell","Hydathode","PP2,XP3","u.a")
  param <- domcell_types(expression_data)
    
  #Define a matrix to store the rows and column counts for each cluster type assay.
  pl <- length(param)
  count_mat <- matrix(pl,ncol = 2)
  exp_found <- c(length=pl)
  
  #Test formula
  for(i in 1:pl)
  {
    #Find cell type.
    cell_type <- param[i]
    atlasRes <- searchAtlasExperiments(properties = cell_type, species = y)
    
    
    if(is.null(atlasRes))
    {
      exp_found[i] <- 0
     # sig <- 0
     # return(sig)
      next
    }
    else
      exp_found[i] <- 1
    
    
    #Get the data.
    atlasData1 <- getAtlasData(atlasRes$Accession[1])
    rnaseqExps1 <- getAtlasData( 
      atlasRes$Accession[1][ 
        grep( 
          "rna-seq", 
          atlasRes$Type, 
          ignore.case = TRUE 
        ) 
      ] 
    )
    
    #Access data.
    allExps1 <- atlasData1[[atlasRes$Accession[1]]]
    sumarrexp1 <- allExps1$rnaseq
    
    assaydata1 <- assays(sumarrexp1)$counts
    head(assaydata1)
    
    rname_assay1 <- rownames(assaydata1)
    colData(sumarrexp1)
    metadata(sumarrexp1)
    allExps <- getAtlasExperiment(atlasRes$Accession[1])
    allExps
    
    dim(assaydata1)
    assay_compute_mat <- as.data.frame(assaydata1)
    
    count_mat[i,1] <- nrow(assay_compute_mat)
    count_mat[i,2] <- ncol(assay_compute_mat)
    
    write.csv(assay_compute_mat,paste("/Domain_Knowledge/Experiment Acession, Cell No.:",i))
  }
  
  #Find the max gene, min gene across all cluster types
  head(expression_data)
    cluster_experiments_details <- c(length=cluster_assignments)
  lp <- length(param)+1
  seuratobj_data <- seuratobj_data_copy
  expression_data <- seuratobj_data[['RNA']]$counts
  rnames <- rownames(expression_data)
  rname_count <- 0
  exist_count <- c(length=length(param))
  nr_exp <- nrow(expression_data)
  nc_exp <- ncol(expression_data)
  geom_mean_clustermat <- matrix(nrow = length(param),ncol = nr_exp)
  count_gene_presence <- matrix(nrow=nr_exp,ncol = length(param))
  pl <- length(param)
  count_mat <- matrix(pl,ncol = 2)
  exp_found <- c(length=pl)
  
  for(i in 1:pl)
  {
    exists <- 1
    cell_type <- param[i]
    
    atlasRes <- searchAtlasExperiments(properties = cell_type, species = "Arabidopsis Thaliana")
    if(is.null(atlasRes))
    {
      geom_mean_clustermat[i,] <- NA
      next
    }  
    
  #  print(paste("Cell Type:",cell_type))
    #cluster_experiments_details[1] <- atlasRes
    
    atlasData1 <- getAtlasData(atlasRes$Accession[1])
    rnaseqExps1 <- getAtlasData( 
      atlasRes$Accession[1][ 
        grep( 
          "rna-seq", 
          atlasRes$Type, 
          ignore.case = TRUE 
        ) 
      ] 
    )
    allExps1 <- atlasData1[[atlasRes$Accession[1]]]
    sumarrexp1 <- allExps1$rnaseq
    
    assaydata1 <- assays(sumarrexp1)$counts
    head(assaydata1)
    
    rname_assay1 <- rownames(assaydata1)
    colData(sumarrexp1)
    metadata(sumarrexp1)
    allExps <- getAtlasExperiment(atlasRes$Accession[1])
    allExps
    
    dim(assaydata1)
    assay_compute_mat <- as.data.frame(assaydata1)
    nr_assay <- nrow(assay_compute_mat)
    
    rowavg <- c()
    for(k in 1:nr_assay)
    {
      cur_row <- assay_compute_mat[k,]
      x <- exp(mean(log(cur_row[cur_row>0])))
      rowavg[k]  <- x
    }
    assay_compute_mat$geom_mean <- rowavg
    #Check how many genes from the seuratobject are present in assaydata1
    #Retreive unnormalized data
    rnames_assay <- rownames(assay_compute_mat)
    fine <- 1
    if(all(rnames_assay==rnames)){  
      print('Perfect match in same order')
    }else if(!all(rnames_assay==rnames) && all(sort(rnames_assay)==sort(rnames))){  
      assay_compute_mat <- assay_compute_mat[order(rnames),]
    }else
      fine <- 0
    
    match_out <- match(rnames_assay,rnames)
    match_positive <- which(!(is.na(match_out)))
    exist_count[i] <- sum(!(is.na(match_out)))
    count_gene_presence[match_positive,i] <- 'Y'
    
    nan_v <- which(is.nan(assay_compute_mat$geom_mean))
    #av <- rowMeans(assay_compute_mat)
    assay_compute_mat$geom_mean[nan_v] <- 0
    
    max_geom <- max(assay_compute_mat$geom_mean)
    min_geom <- min(assay_compute_mat$geom_mean)
    
    if(!(fine==0)){
      geom_mean_clustermat[i,] <- assay_compute_mat$geom_mean
    }else{
      not_nan <- which(!is.na(match_out))
      geom_mean_clustermat[i,not_nan] <- assay_compute_mat$geom_mean[not_nan] 
    }
    
    #Write all relevant supplementary data to a file.
    write.csv(assay_compute_mat,paste(i,paste("/Dom_knowledge/Assay_compute_Matrix.csv",i)))
    write.csv(colData(sumarrexp1),paste("/Dom_knowledge/ColNames_Summary.csv",i))
  }
  names(exist_count) <- param
  
  #write.csv(exist_count,"Count of Genes Found in the 19 Cluster Cell Types.")
  #write.csv(count_gene_presence,paste(i,"Count_Gene_AssayMatrix.csv"))
  sigmoid_fuzz <- c(length=nr_assay)
  top_quan <- quantile(assay_compute_mat$geom_mean)["75%"]
  top_quan <- c()
  
  #Find row quantiles
  library(matrixStats)
  top_quan <- matrix(nrow = length(param),ncol = nr_exp)
  trans_geom_mean_clustermat <- t(geom_mean_clustermat)
  top_quan <- rowQuantiles(trans_geom_mean_clustermat,na.rm = TRUE)
  top_25_quan <- top_quan[,2]
  sigmoid_membership_mat <- matrix(nrow = nr_exp,ncol = length(param))
  
  #Calculate membership of genes in the 19 clusters
  for(i in 1:nr_exp)
  {
    for(j in 1:length(param))
    {
      x <- trans_geom_mean_clustermat[i,j]
      if(is.na(x))
      {
        sigmoid_membership_mat[i,j] <- NA
      }
      else
      {
        sigmoid_membership_mat[i,j] <- 1/(1+exp(-(x-top_25_quan[i])))
      }
    }
  }
  write.csv(cbind(rnames,sigmoid_fuzz),"/Dom_knowledge/Sigmoid_Count_Values.csv")
  rownames(geom_mean_clustermat) <- param
  colnames(geom_mean_clustermat) <- rnames
  
  rownames(trans_geom_mean_clustermat) <- rnames
  colnames(trans_geom_mean_clustermat) <- param
  
  rownames(sigmoid_membership_mat) <- rnames
  colnames(sigmoid_membership_mat) <- param
  #Write sigmoid membership function 
  
  #Write the trans_geom_matrix and sigmoid function into csv files
  write.csv(trans_geom_mean_clustermat,"/DomKnowledge/TransGeom_GSE161332.csv")
  write.csv(sigmoid_membership_mat,"/DomKnowledge/Sigmoid_Membership_Clusters.csv")
  return(sigmoid_membership_mat)
}


#' Celsius conversion
#'
#' Find the cell types from any expression matrix
#' @param X The expression matrix
#' @return The compressed matrix
#' @examples 
#' temp1 <- compRes(22);

#' @export

normRes <- function(x)
{
  
  expression_data <- x
  clus_type <- domcell_types(x)
  names_gene <- rownames(expression_data)
  lclus <- length(clus_type)
  all_exist <- c(lclus)
  max_clusgene <- c(lclus)
  min_clusgene <- c(lclus)
  
  for(i in 1:19)
  {
    if(!file.exists(paste("Experiment Acession, Cell No.:",i)))
    {
      max_clusgene[i] <- "not found"
      min_clusgene[i] <- "not found"
      next
    }
    mat <- read.csv(paste("Experiment Acession, Cell No.:",i))
    rownames(mat) <- mat[,"X"]
    mat <- mat[,-1]
    rmat <- rownames(mat) 
    if(all(names_gene %in% rmat))
      all_exist[1] <- 1
    mat_subset <- mat[names_gene, ] 
    
    mat_subset <- mat[rmat %in% names_gene, ]  
    rmed <- rowMedians(as.matrix(mat_subset))
    
    max_i <- which.max(rmed)
    min_i <- which.min(rmed)
    
    max_clusgene[i] <- names(max_i)
    min_clusgene[i] <- names(min_i)
    
  }
  
  #Print the names of max genes 
  print(max_clusgene)
  print(min_clusgene)
  max_g <- names(sort(table(max_clusgene),decreasing=TRUE)[1])
  min_g <- names(sort(table(min_clusgene),decreasing = TRUE)[1])
  
  #Find the maximum number of cells across all clusters
  max <- -1
  for(i in 1:19)
  {
    s <- sum(cluster_assignments==i)
    if(max<s)
      max <- s
    
  }
  
  #Take a subset of expression data and try out
  #Create a 3D matrix to store genes vs. cells w.r.t. clusters
  genecClus <- array(dim = c(length(names_gene),max,19))
  
  #Normalize by max across all clusters
  ind_m <- which(names_gene==max_g)
  lim <- ind_m+5
  exp_subs <- expression_data[ind_m:lim,]
  exp_subs_norm <- exp_subs
  #Scaling operation
  for (i in 2:nrow(exp_subs)) {
    for(j in 1:ncol(exp_subs))
    {
      if(exp_subs[i,j] < exp_subs[1,j])
      {
        
        scaling_fac <- exp_subs[i,j]/(1/max(exp_subs[,j])*abs(exp_subs[1,j]-exp_subs[i,j])+1)
        exp_subs_norm[i,j] <- exp_subs[i,j]*scaling_fac
      }
      else if(exp_subs[i,j] > exp_subs[1,j])
      {
        scaling_fac <- exp_subs[i,j]*(1/max(exp_subs[,j])*abs(exp_subs[1,j]-exp_subs[i,j])+1)
        exp_subs_norm[i,j] <- exp_subs[i,j]*scaling_fac
      }
    }
  }
  
  ratio <- c()
  for(i in 1:nrow(expression_mat))
  {
    sp <- sum(expression_mat[i,]==0)
    len <- length(expression_mat[i,])
    
    ratio[i] <- sp/len
  }
  av_rat <- mean(ratio)
  
  
  ##Normalize by max per cluster
  genecClus_norm <- array(dim = c(5,max,19))
  for(i in 1:19)
  {
    clus_no <- which(cluster_assignments==i)
    
    ngene <- max_clusgene[i]
    if(ngene=="not found")
      next
    ind_gene <- which(names_gene==ngene)
    maxclus_row <- expression_data[ind_gene,clus_no]
    
    genecClus[1,1:length(clus_no),] <- maxclus_row
    low <- ind_gene+1
    high <- ind_gene+4
    x <- 2
    for(j in low:high)
    {
      genecClus[x,1:length(clus_no),i] <- expression_data[j,clus_no]
      x <- x + 1
    }
  }
  
  exclude_ind <- which(max_clusgene=="not found")
  genecClus_norm <- array(dim = c(5,max,19))
  
  for (i in 1:19) {
    if(i %in% exclude_ind)
      next
    ncol <- max
    genecClus_norm[1,,i] <- genecClus[1,,i]
    for(j in 2:5)
    {
      na_col <- which(is.na(colSums(genecClus[,,i])))[1]
      
      if(!(is.na(na_col)))
        ncol <- na_col-1
      
      for(k in 1:ncol)
      {
        mcol <- max(genecClus[,k,i])
        if(genecClus[j,k,i] < genecClus[1,k,i])
        {
          scaling_fac <- genecClus[j,k,i]/((1/mcol)*abs(genecClus[1,k,i]-genecClus[j,k,i])+1)
          genecClus_norm[j,k,i] <- genecClus[j,k,i]*scaling_fac
        }
        else if(genecClus[j,k,i] > genecClus[1,k,i])
        {
          scaling_fac <- genecClus[j,k,i]*((1/mcol)*abs(genecClus[1,k,i]-genecClus[j,k,i])+1)
          genecClus_norm[j,k,i] <- genecClus[j,k,i]*scaling_fac
        }
      }
    }
  }
  
  #Now combine all the clusters into one matrix
  genecClusnorm <- matrix(nrow = 5,ncol = ncol(expression_data))
  for(k in 1:19)
  {
    if(k %in% exclude_ind)
      next
    na_col <- which(is.na(colSums(genecClus[,,k])))[1]
    
    if(!(is.na(na_col)))
      ncol <- na_col-1
    
    start_col <- 1
    for(i in 1:5)
    {
      genecClusnorm[i,start_col:ncol] <- genecClus_norm[i,1:ncol,k]
    } 
    start_col <- ncol+1
  }
  
  #Write all the csv files
  library(xlsx)
  #Across All Clusters
  write.xlsx2(exp_subs, "AllClus_Norm.xlsx", sheetName = "BeforeNorm",
              col.names = TRUE, row.names = TRUE, append = TRUE)
  
  write.xlsx2(exp_subs_norm, "AllClus_Norm.xlsx", sheetName = "AfterAllClusNorm",
              col.names = TRUE, row.names = TRUE, append = TRUE)
  
  #Per Cluster
  write.xlsx(genecClus_norm, "AllClus_Norm.xlsx", sheetName = "AfterPerClusNorm",
             col.names = TRUE, row.names = TRUE, append = TRUE)
  
  exp_subsSD <- apply(exp_subs, 1, sd)
  exp_subsMean <- rowMeans(exp_subs)
  
  exp_subs_normSD <- apply(exp_subs_norm, 1, sd)
  exp_subs_normMean <- rowMeans(exp_subs_norm)
  
  exp_genecClus_SD <- apply(genecClusnorm, 1, sd)
  exp_genecClus_Mean <- rowMeans(genecClusnorm)
  
  library(vsn)
  meanSdPlot(as.matrix(exp_subs))
  meanSdPlot(as.matrix(exp_subs_norm))
  meanSdPlot(genecClusnorm)
  
  t_zscore <- read.csv("C:\\Users\\user\\Desktop\\Research\\SCRNA-NORMALISATION\\Test_zscore.csv")
  t_zscore
  
  tzscore_SD <- apply(t_zscore, 1, sd)
  tzscore_Mean <- rowMeans(t_zscore)
  
  meanSdPlot(as.matrix(t_zscore))
  
  
  ind_m <- which(names_gene==max_g)
  lim <- ind_m+5
  exp_subs <- expression_data[ind_m:lim,]
  exp_subs_norm <- exp_subs
  #Scaling operation
  for (i in 2:nrow(exp_subs)) {
    for(j in 1:ncol(exp_subs))
    {
      if(exp_subs[i,j] < exp_subs[1,j])
      {
        
        scaling_fac <- exp_subs[i,j]/(1/max(exp_subs[,j])*abs(exp_subs[1,j]-exp_subs[i,j])+1)
        exp_subs_norm[i,j] <- exp_subs[i,j]*scaling_fac
      }
      else if(exp_subs[i,j] > exp_subs[1,j])
      {
        scaling_fac <- exp_subs[i,j]*(1/max(exp_subs[,j])*abs(exp_subs[1,j]-exp_subs[i,j])+1)
        exp_subs_norm[i,j] <- exp_subs[i,j]*scaling_fac
      }
    }
    
  }
  
}
  
  
  

  
  
  
  
  
  
  
  
  
