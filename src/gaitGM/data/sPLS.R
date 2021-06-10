######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
#
# DESCRIPTION: Perform a sparse PLS over subsets of gene Expression and metabolite data.
#
# VERSION: 1.0
#######################################################################################

sPLS <- function(geneData, metData, keepX){
  "
  This function perform an sparse PLS (sPLS) using mixOmics spls function. Genes are the
  explanatory variable (X), while metabolites are the response variables (Y).

  Arguments:
    :param geneData: Gene Expression Wide Dataset.
    :type geneData: Matrix

    :param metData: Metabolite Wide Dataset.
    :type metData: Matrix

    :param keepX: Number of genes to keep in each component for creating the model.
    :type keepX: integer

  Returns:
    :return SPLS: spls object from the function spls().
    :type SPLS: spls object
  "

  # Libraries
  library("mixOmics")
  library("RColorBrewer")

  # Force to numbers
  indx <- sapply(geneData, is.factor)
  geneData[indx] <- lapply(geneData[indx], function(x) as.numeric(as.character(x)))
  indx <- sapply(metData, is.factor)
  metData[indx] <- lapply(metData[indx], function(x) as.numeric(as.character(x)))

  # Transpose matrices (Genes in columns, Subjects in rows)
  geneData_T <- t(geneData)
  metData_T <- t(metData)

  # sPLS
  SPLS <- spls(geneData_T, metData_T, ncomp = 2, keepX = c(as.integer(keepX),as.integer(keepX)), scale = TRUE, 
               mode = "classic", near.zero.var = TRUE)

  return(SPLS)
}

plotInPdf <- function(splsObjects, figurePath, multipleNames){
  "
  This function uses the cim() function from mixOmics to plot a heatmap with the sPLS results.

  Arguments:
    :param splsObjects: List with all the SPLS objects from sPLS function.
    :type splsObects: List

    :param figurePath: Full path for the sPLS output figure.
    :type figurePath: string

    :param multipleNames: List containing the different metabolite types (e.g. sphingomyelin).
    :type multipleNames: List
  "

  pdf(file=figurePath, height = 12, width = 10)

  i <- 1
  for (SPLS in splsObjects){
    cim(SPLS,
        row.cex = 0.75, 
        col.cex = 0.75,
        transpose = FALSE,
        mapping = "XY",
        margins = c(8,12.5),
        title = multipleNames[i],
        legend = list(title = "EXPLAINED VARIANCE",
                      legend = list(paste("Comp 1: ", round(SPLS$explained_variance$Y[1]*100, 0), '%'), 
                                    paste("Comp 2: ", round(SPLS$explained_variance$Y[2]*100, 0), '%')
                                    ),
                      col = 'grey', cex = 0.8),
        dist.method = c("correlation", "correlation")
        )
    print(paste(multipleNames[i], 'done.', sep = " "))
    i <- i+1
  }
  dev.off()
}

corrMat <- function(splsObjects, multipleNames, threshold){
  "
  This function creates again the correlation matrix between genes and metabolites and creates a sif-like
  table with the most important correlations.

  Arguments:
    :param splsObjects: List with all the SPLS objects from sPLS function.
    :type splsObects: List

    :param multipleNames: List containing the different metabolite types (e.g. sphingomyelin).
    :type multipleNames: List

    :param threshold: PValue threshold. Only values under this value will be shown in the results table.
    :type threshold: float

  Returns:
    :return correlations: sif-like file with this information: Metabolite Gene Correlation Subset
    :rtype correlations: Matrix
  "

  correlations <- as.data.frame(t(c("Metabolite", "Gene", "Correlation", "Subset")))
  colnames(correlations) <- correlations[1,]
  correlations <- correlations[-1,]
  i <- 1

  for (SPLS in splsObjects){
    comp=1:SPLS$ncomp

    keep.X = apply(abs(SPLS$loadings$X[,comp, drop = FALSE]), 1, sum) > 0
    keep.Y = apply(abs(SPLS$loadings$Y[,comp, drop = FALSE]), 1, sum) > 0

    cord.X = cor(SPLS$X[, keep.X, drop = FALSE], SPLS$variates$X[, comp], use = "pairwise")
    cord.Y = cor(SPLS$Y[, keep.Y, drop = FALSE], SPLS$variates$X[, comp], use = "pairwise")

    XY.SPLS = as.matrix(cord.X %*% t(cord.Y))

    for(row in 1:nrow(XY.SPLS)) {
      for(col in 1:ncol(XY.SPLS)) {
        correlation <- as.numeric(round(XY.SPLS[row,col], digits = 3))
        if (abs(correlation) > threshold){
          metabolite <- colnames(XY.SPLS)[col]
          gene <- rownames(XY.SPLS)[row]
          newLine <- data.frame(metabolite, gene, as.numeric(correlation), multipleNames[i], stringsAsFactors = FALSE)
          names(newLine) <- names(correlations)
          correlations <- rbind(correlations, newLine)
        }
      }
    }
    i <- i+1
  }
  if(nrow(correlations) == 0){
    correlations <- paste("There are no stronger correlations over", threshold, "- Try a smaller threshold.")
  }else{
    colnames(correlations) <- c("Metabolite", "Gene", "Correlation", "Subset")
    correlations <- correlations[order(correlations[4], abs(correlations[3]), decreasing = TRUE), ]
  }

  return(correlations)
}
