##### Simple script to verify PANA on Allison data

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





# Scripts and variables
  source("PCA2GO.2.R")  # function
  source("PCA-GENES.R") # function
  PCA2GO.fun = "PCA2GO"  # change to PCA2GO.2 when  B2GScore was not run
  sel = "single%"  # this is the selection criteria for a metagene, means it is based on a single PC (this can be changed)
  cutoff =  0.3 # the single PC should explain at least 30% of the varianze to be captured (this can be changed)
  cutoff = 2
  sel = "rel.abs"
  
##### Get data
input_data <- read.delim ("wide_genes.tsv", as.is = TRUE, row.names = 1)  # expression data
head (input_data)  # we see that gene names are in the rownames of the matrix

genes2pathway <- read.delim ("pathway_gene.tsv", as.is = TRUE)
head (genes2pathway)  # this file contains multiple columns, we only need geneID and pathway ID.
annotation <- data.frame(gene = genes2pathway$UniqueID, path= genes2pathway$Pathway_Name)

# Now we need the same way to naming genes in ghe input.data and the annotation file. We need to add "gene_" to the gene number
#rownames(input_data) <- paste("gene_", rownames(input_data), sep = "") # needed to match the genes
#head(input_data)  # Fixed!

# Running PANA analysis -----------------------------------------------------------------------------------------
expression_GO <- get(PCA2GO.fun)(input_data, annotation, var_cutoff = cutoff, fac_sel =  sel)
metagenes <- expression_GO$X.sel
head(metagenes)



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

#  Visualizing metagenes
library(RColorBrewer)
library("d3heatmap")
heatmap(metagenes, Colv = NA, keep.dendro = FALSE, col = rev(brewer.pal(9,"RdBu")))
scaled.data<- t(scale(t(metagenes)))
d3heatmap(t(scale(t(metagenes))), colors = rev(brewer.pal(9,"RdBu")),Colv = FALSE,
          dendrogram = 'none', Rowv = FALSE, cexRow = 0.9, cexCol = 0.4, 
          yaxis_width=300)
