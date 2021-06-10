######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
#
# DESCRIPTION: Perform a correlation analysis of two datasets.
#
# VERSION: 1.0
#######################################################################################

#################
# Main Function #
#################

corr_main_func <- function(x, y, meth, thres, corrMatPath, outputPath, figurePath){

  "
  This function performs a correlation analysis between 2 matrices (genes and metabolites).
  For obtaining the PValue of the correlation, it simulates normal data for both original datasets
  1000 times.

  Arguments:
    :param x: Gene Wide Dataset with samples in columns and genes in rows.
    :type x: matrix

    :param x: Metabolite Wide Dataset with samples in columns and metabolites in rows.
    :type x: matrix

    :param meth: Methodology for the correlation function. One of 'pearson', 'spearman' or 'kendall'.
    :type meth: string

    :param thres: Threshold to cut the correlations for the output table.
    :type thres: float

    :param corrMatPath: Full path for the Correlation Matrix output
    :type corrMatPath: string

    :param outputPath: Full path for the Output table
    :type outputPath: string

    :param figurePath: Full Path for the Network-like output figure
    :type figurePath: string
  "

  # Transpose matrices
  t_x <- t(x)
  t_y <- t(y)
  # Check NearZeroVar
  new_x <- checkZeroVar(t_x)
  new_y <- checkZeroVar(t_y)
  ## Estandarize # Finally Removed
  # new_x <- scale(new_filter_x)
  # new_y <- scale(new_filter_y)
  # Original correlation
  corr_orig <- cor(new_x, new_y, method = meth)
  # Number of iterations
  N <- 1000
  # Create a pvalue matrix
  pval_mat_counts <- corr_orig
  pval_mat_counts[,] <- 0
  # Total Number of correlations
  values <- nrow(corr_orig)*ncol(corr_orig)
  # Vector of means
  means_x <- apply(new_x, 2, mean)
  means_y <- apply(new_y, 2, mean)
  # Vector of sds (1s if scaled)
  sds_x <- apply(new_x, 2, sd)
  sds_y <- apply(new_y, 2, sd)
  for (i in 1:N){
    # New datasets
    x_new <- create_newDataset(new_x, means_x, sds_x)
    y_new <- create_newDataset(new_y, means_y, sds_y)
    corr_new <- cor(x_new, y_new, method = meth)
    for(j in 1:values){
      new_corr_value <- corr_new[j]
      orig_corr_value <- corr_orig[j]
      if (abs(new_corr_value) > abs(orig_corr_value)){
        pval_mat_counts[j] <- pval_mat_counts[j] + 1
      }
    }
  }
  # Final pvalue matrix
  pval_mat <- pval_mat_counts/N
  finalTable <- return_result(corr_orig, pval_mat, thres)

  # Write outputs
  plot_relations_network(finalTable, figurePath)
  write.table(corr_orig, corrMatPath, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
  write.table(finalTable, outputPath, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#######################
# Auxiliary Functions #
#######################

# Check Near Zero Variance
checkZeroVar <- function(matrix){

  "
  This function uses the nearZeroVar function from caret to remove the variables (genes or metabolites)
  with near zero variance.

  Arguments:
    :param matrix: Input matrix to analyze the nzv.
    :type matrix: matrix

  Returns:
    :return newMatrix: Matrix without near zero variance variables.
    :rtype newMatrix: matrix
  "
  library(caret)

  nzv_cols <- nearZeroVar(matrix)
  if (length(nzv_cols) != 0){
    newMatrix <- matrix[,-nzv_cols]
  }
  else{
    newMatrix <- matrix
  }
  return(newMatrix)
}

create_newDataset <- function(dataset, mean_list, sd_list){

  "
  This function simulates a new dataset with random values following a normal with different mean for
  each feature (gene or metabolite) and standard deviation of 1.

  Arguments:
    :param dataset: original matrix to simulate values.
    :type dataset: Matrix

    :param mean_list: list with the means for each feature (gene or metabolite)
    :type mean_list: list

    :param sd_list: list with the standard deviations for each feature (gene or metabolite)
    :type sd_list: list

  Returns:
    :return newDataset: Simulated dataset with normal distribution for each feature (gene or metabolite)
    :rtype bewDataset: Matrix
  "

  newDataset <- mapply(rnorm, nrow(dataset), mean_list, sd_list)

  return(newDataset)
}

return_result <- function(corr_mat, pval_mat, thres){

  "
  This function returns the All vs All correlation result in format: Gene Metabolite Correlation Pvalue.
  One value is only written if the PValue is smaller than the stablished threshold.

  Arguments:
    :param corr_mat: Correlation matrix between genes and metabolites.
    :type corr_mat: Matrix

    :param pval_mat: PValue matrix for the correlations between genes and metabolites.
    :type pval_mat: Matrix

    :param thres: PValue threshold
    :type thres: float

  Returns:
    :return finalResult: Output table with this information: Gene Metabolite Correlation Pvalue.
    :rtype finalResult: sif-like file
  "

  values <- nrow(corr_mat)*ncol(corr_mat)
  finalResult <- data.frame("gene", "metabolite", "correlation", "p-value", stringsAsFactors = FALSE, check.names = FALSE)
  finalResult <- finalResult[-1,]
  for (i in 1:values){
    if (pval_mat[i] <= thres){
      k <- arrayInd(i, dim(pval_mat))
      gene <- rownames(pval_mat)[k[,1]]
      met <- colnames(pval_mat)[k[,2]]
      newLine <- data.frame(gene, met, as.numeric(corr_mat[i]), as.numeric(pval_mat[i]), stringsAsFactors = FALSE)
      names(newLine) <- names(finalResult)
      finalResult <- rbind(finalResult, newLine)
    }
  }
  finalResult <- finalResult[order(abs(finalResult[3]), decreasing = TRUE),]

  return(finalResult)
}

plot_relations_network <- function(relations.sif, figurePath){

  "
  This function uses igraph for plotting the correlations in a network-like file.

  Arguments:
    :param relations.sif: finalResult table from return_result function.
    :type relations.sif: sif-like file

    :param figurePath: Full path for saving the pdf output.
    :type figurePath: string
  "

  # Libraries
  library(igraph)

  # Only plot first 500 strong correlations
  numRelations <- nrow(relations.sif)
  if (numRelations>500){
    relations.sif <- relations.sif[c(1:500),]
    numRelations <- 500
  }

  graph <- graph.data.frame(relations.sif, directed = FALSE)
  V(graph)[relations.sif[,1]]$color <- "royalblue1"
  V(graph)[relations.sif[,2]]$color <- "forestgreen"
  E(graph)$color <- ifelse(relations.sif[,3] > 0,'firebrick2','deepskyblue')
  pdf(file=figurePath, height = 11, width = 8.5)

  plot.igraph(graph,
              layout=layout.fruchterman.reingold,
              vertex.size=4,
              vertex.frame.color="gray",
              vertex.label.color="black",
              vertex.label.cex=0.3)
  legend("bottomright",
         legend = c("Genes", "Metabolites", "Positive Correlation", "Negative Correlation"),
         col = c("royalblue1","forestgreen", 'firebrick2', 'deepskyblue'),
         pch = c(19,19,95,95),
         bty = "n",
         pt.cex = 1,
         cex = 0.7,
         text.col = "black",
         horiz = F)
  title(main = paste(as.character(numRelations), "Strongest correlations\nBetween Genes and Metabolites", sep = " "))

  dev.off()
}
