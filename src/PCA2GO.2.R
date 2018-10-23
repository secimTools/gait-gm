##############################################################
# This functions performs PCA and selects scores for genes
# associated by belonging to a GO, given a selection of GOs
# It prodices a matrix X.sel with columns the number of
# conditions and rows the number of GO_genes
#
# Input data
# X: expression data matrix
# selection: GO selection object obtained by select.GO
# fac_sel = criterion to select components, con be:
#   %accum = percentage of accumulated variability
#   single% = percentage of variability of that PC
#   abs.val = absolute value of the variabily of that PC
#   rel.abs = fold variability of tot.var/rank(X)
# var_cutoff = variability cut off value
#
# Ana Conesa aconesa@cipf.es 8 September 2007
############################################################

PCA2GO <- function (X, annotation, var_cutoff, fac_sel)
{
var_cutoff <- as.numeric(var_cutoff)
fac_sel <- match.arg(fac_sel, c("accum", "single", "abs.val", "rel.abs"))
# initialize variables
X.sel <- NULL 
n.go <- NULL
total.genes <- NULL
variab <- vector (mode="numeric")
tot.variab <- vector (mode="numeric")
eigen.val <- NULL
go.sel <- unique(annotation[,2])
n.ge <- NULL
X.loadings <- vector(mode = "list", length = 0)
# PCAs for all GOs loop
   for (i in 1: length(go.sel)) {
    gene.sel <- annotation[annotation[,2] == go.sel[i],1]
        if (length(gene.sel) > 1) {
      total.genes <- unique(c(total.genes, gene.sel))
      gene.sel <- is.element(rownames(X), gene.sel)
      if (length(which(gene.sel)) > 1) {
         gene.sel <- X[gene.sel,]
         if (any(is.na(gene.sel))) { gene.sel <- pca4NA(gene.sel)}
        pca.sel <- PCA.GENES(t(gene.sel))  # pca
        eigen <- pca.sel$eigen$values
        tot.var <- sum(eigen)
        eigen.val <- c(eigen.val, tot.var)
        rank <- length(which(eigen > 1e-16))
        level <- 1
        # num fac
          if (fac_sel == "accum") {
            fac <- max(length(which(pca.sel$var.exp[,2] <= var_cutoff / sqrt(level))),1)
        } else if (fac_sel == "single"){
              fac <- length(which(pca.sel$var.exp[,1] >= (var_cutoff / sqrt(level))))
        } else if (fac_sel == "rel.abs"){
            mean.expl.var <- tot.var/ nrow(gene.sel)
            fac <- length(which(eigen >= (mean.expl.var*var_cutoff / sqrt(level))))
        } else if (fac_sel == "abs.val"){
            abs.val.bycomp <- mean(apply(pca.sel$Xoff,2,var)) * nrow(gene.sel)
            fac <- length(which(eigen >= abs.val.bycomp * var_cutoff / sqrt(level)))
        }
        tot.variab <- c(tot.variab, pca.sel$var.exp[,1])
        #variab <- c(variab, pca.sel$var.exp[1:fac,1])
        if (fac > 0 ) { # num fac
                variab <- c(variab, pca.sel$var.exp[1:fac,1])
                n.ge <- c(n.ge, nrow(gene.sel))
              n.go <- c(n.go, fac)
              data.h <- as.matrix(pca.sel$scores[,1:fac])
                loads <-  as.matrix(pca.sel$loadings[, 1:fac])
        colnames(data.h) <- paste(go.sel[i], c(1:fac), sep="_")
              X.sel <- cbind(X.sel,data.h) # attach to results matrix
                for (u in 1:fac)  { X.loadings[[length(X.loadings)+1]] <- loads[,u] }
                   }
            }
    }
   }
# Rearrange result
rownames(X.sel) <- colnames(X)
names(X.loadings) <- colnames(X.sel)
X.sel <- t(X.sel)
# Paco edit (Only X.sel needed)
expression_GO <- list("X.sel" = X.sel, "X.loadings" = X.loadings,"n.ge" = n.ge,"n.go" = n.go, "go.sel" = go.sel, "total.genes" = total.genes, "var_cutoff" = var_cutoff, "variab" = variab, "tot.variab" = tot.variab, "eigen.val" = eigen.val)
metagenes <- as.data.frame(expression_GO$X.sel)

# ADJUST GENE DIRECTION -----------------------------------------------------------------------------------------
##for each metagene
for (i in 1:length(row.names(metagenes)) ) {
  cur_pathway_id <- unlist(lapply(strsplit(row.names(metagenes)[i], "_"), function(x) x[1]))
  #loadings indicates the contribution of each gene to PC1
  #The loading sign is not arbitrary.
  #Positive loading indicates positive correlation of gene expression with the scores while negative loading indicates negative correlation.
  gene_loadings <- unlist(expression_GO$X.loadings[i])
  #select genes that contribute most in a given component as
  # abs(loading del gen)/sum(abs(loadings de todos los genes)
  nGenes <- length(gene_loadings)
  loadings_sum <- sum(abs(gene_loadings))
  has_positives<-0
  has_negatives<-0
  selected <- c()
  for(j in gene_loadings ){
    #If this value is greater than 1/total_genes, the gene is selected because it has a greater contribution
    #than the value of contribution if all the same genes contribute together.
    if(abs(j)/loadings_sum > 1/nGenes){
      selected <- c(selected, j)
      has_positives<- has_positives + ifelse(j > 0, 1, 0)
      has_negatives<- has_negatives + ifelse(j < 0, 1, 0)
    }
  }

  ## Change the direction for the metagene
  ##If most or all of the genes have a positive loading then
  if(has_positives > has_negatives ){
    #leave metagene as it is
    ##If most or all of the genes have a negative loading then invert metagene
  }else if(has_negatives > has_positives){
    metagenes[i,] <-metagenes[i,] * -1
    ##If same number of negative and positive loadings then resolve
  }else if(has_negatives > 0 && has_positives > 0){
    has_positives <- sum(selected[selected>0])
    has_negatives <- abs(sum(selected[selected<0]))
    if(has_negatives > has_positives){
      ##If negative loadings genes are bigger then invert metagene
      metagenes[i,] <-metagenes[i,] * -1
    }
  }
}

# Output
PANAOutput <- preparePANAOutput(X, X.loadings)
# ri2py is not taking rownames ??? - Solved
metagenes <- cbind(rownames(X.sel), metagenes)
colnames(metagenes)[1] <- "metagene_name"
PANAOutput <- cbind(rownames(X), PANAOutput)
colnames(PANAOutput)[1] <- "KEGG_ID"
output <- list(metagenes, PANAOutput)
return(output)
}

PCA.GENES<-function(X)
{
  #PCA.GENES is very useful to obtain principal components to a matrix that has more variables than individuals. 
  #R can not apply princomp is such case and when there are a lot of variables eigen(t(X)%*%X) can not be computed.

  #X is a matrix that has on columns the genes considered as variables in the PCA analysis.
  #First we center the matrix by columns (Xoff) and then we obtain the eigenvalues and the eigenvectors of the matrix Xoff%*%t(Xoff) and we #use the equivalences between the loadings and scores to obtain the solution
  #Llamo scores1 y loadings1 a lo que busco y scores2 y loadings2 a los scores y loadings de la traspuesta

  X <- as.matrix(X)
  n<-ncol(X)
  p<-nrow(X)
  offset<-apply(X,2,mean)
  Xoff<-X-(cbind(matrix(1,p,1))%*%rbind(offset))

  #eigen command sorts the eigenvalues in decreasing orden.

  eigen<-eigen(Xoff%*%t(Xoff)/(p-1))
  var<-cbind(eigen$values/sum(eigen$values),cumsum(eigen$values/sum(eigen$values)))

  loadings2<-eigen$vectors
  scores2<-t(Xoff)%*%loadings2

  normas2<-sqrt(apply(scores2^2,2,sum))

  scores1<-loadings2%*%diag(normas2)
  loadings1<-scores2%*%diag(1/normas2)

  output<-list(eigen,var,scores1,loadings1, Xoff)
  names(output)<-c("eigen","var.exp","scores","loadings", "Xoff")
  output
}

pca4NA <- function (X, fac = min(dim(X)), conver = 1e-07, max.iter = 1000) 
{
  if (any(is.na(X))) { 
    X <- as.matrix(X)
    print("Missing values are taken care of by PCA")
    NA.position <- which(is.na(X)) # finds NAs of X
    print(paste(length(NA.position), "missing values"))
    NA.values <- rnorm(length(NA.position), mean(X, na.rm = T), 
                       sd(X, na.rm = T)) # creates random values
    X[NA.position] <- NA.values # put random values in X
    SST <- 0
    for (it in 1:max.iter) { # start loop
      SST.old <- SST
      pca <- PCA.GENES(t(X)) # fits PCA
      T <- as.matrix(pca$scores[,1:fac]) # scores 
      P <- as.matrix(pca$loadings[,1:fac]) # loadings
      Xe <- T %*% t(P) # model
      Xe <- t(Xe)
      X[NA.position] <- Xe[NA.position] # put estimated values in X
      SST <- sum(X^2)
      # Convergence
      if (abs(SST-SST.old)/SST < conver) break
    }
  } else { print("No missing values, unchanged output") }
  return(X) # output
}

preparePANAOutput <- function(expressionMatrix, loadingsList){
  PANAOutput <- data.frame(rownames(expressionMatrix), stringsAsFactors = FALSE)
  rownames(PANAOutput) <- PANAOutput[,1]
  PANAOutput <- PANAOutput[,-1]
  i <- 1
  for (metagene in loadingsList){
    genesList <- as.list(names(metagene))
    PANAOutput[i] <- rep(0, nrow(PANAOutput))
    colnames(PANAOutput)[i] <- names(loadingsList)[i]
    for (gene in genesList){
      PANAOutput[rownames(PANAOutput) == gene, i] <- 1
    }
    i <- i+1
  }
  return(PANAOutput)
}
