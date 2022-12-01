#' @title Perform Gene Expression Significance Test
#'
#' 
#' @param object an ENIGMA object
#' @param nMC the number of bootstrapping times
#' @param p_threshold the threshold of pvalue to determine the significant gene
#' @param refine whether perform refinement
#' @param auto automatically determine whether perform refinement
#' @param filtering perform gene filtering for enigma object. Default: TRUE
#'
#' @return A list object contain two object
#'  call: A binary matrix (row: gene, column: cell type) indicating which gene has high confidence to be accurated estimated
#'  pval: A p-value matrix (row: gene, column: cell type) indicating the confidence
#'
#'
#'
#' @export
GeneSigTest <- function(object,nMC = 1000,p_threshold = 0.05,refine=FALSE,auto=TRUE,filtering=FALSE){
  Bulk = object@bulk
  frac = object@result_cell_proportion
  ref = object@ref
  
  if(auto){
    cor <- cor(ref)
	diag(cor) <- 0
	if(sum(cor>0.85) > 0){
	  refine = FALSE
	}else{
	  refine = TRUE
	}
  }
  
  writeLines("Using Linear Regression To Infer Probability of Expression...")
  p_tab <- NULL
  df <- (nrow(frac)-ncol(frac))
  for(i in 1:nrow(Bulk)){
   exp <- Bulk[i,] 
   rlm.o <- lm(exp ~ as.matrix(frac)-1)
   p <- pt(summary(rlm.o)$coef[,3], df,lower.tail=FALSE)
   p_tab <- rbind(p_tab,p)
  }
  p_tab[is.nan(p_tab)] <- 1
  pvalue.m <- p_tab
  for(i in 1:ncol(pvalue.m)) pvalue.m[,i] <- p.adjust(pvalue.m[,i],method="BH")
  colnames(pvalue.m) <- colnames(frac)
  rownames(pvalue.m) <- rownames(Bulk)
  ### soft max transformation
  if(refine){
  expp <- softmaxP(pvalue.m)
  rownames(expp) <- rownames(pvalue.m) <- rownames(Bulk)
  colnames(expp) <- colnames(pvalue.m) <- colnames(frac)
  
  score <- matrix(NA,nrow = nrow(expp),ncol = ncol(expp))
  rownames(score) <- rownames(expp)
  colnames(score) <- colnames(expp)
  
 
  writeLines("Bootstrapping...")
  for(i in 1:ncol(score)){
   ss <- expp[,i]/rowSums(expp[,-i])
   ss[is.nan(ss)] <- 0
   
   ### using bootstrapping to calculate pvalue
   tab <- list()
   nf <- nrow(p_tab)
   for(n in 1:nMC){
     ss_pe <- expp[,i]/(rowSums(expp[,-i])[sample(1:nf,nf,replace=FALSE)])
	 tab[[n]] <- ss_pe
   }
   tab <- t(matrix(unlist(tab), nrow=length(expp[,1])))
   
   vec <- NULL
   for(j in 1:nf){
      vec <- c(vec,sum(ss[j]<tab[,j]))
   }
   score[,i] <- vec/nMC
  }
  }else{
  score <- pvalue.m
  }
  ####### refine the results and only assign the gene to the most significant cell type
  if(refine) writeLines("Refining...")
  call <- matrix(0,nrow = nrow(score),ncol = ncol(score))
  rownames(call) <- rownames(score)
  colnames(call) <- colnames(score)
  
  for(ct in 1:ncol(score)){
   gene <- rownames(score)[score[,ct]<p_threshold]
   if(refine) order <- apply(pvalue.m[gene,],1,which.min)
   if(refine) gene <- gene[order == ct]
   call[gene,ct] <- 1
  }
  gene.call <- rowSums(call)
  
  ## return the results and filtered enigma object
  if(filtering){
  Bulk = Bulk[gene.call>0,]
  ref = ref[gene.call>0,]
  egm@bulk <- Bulk
  egm@ref <- ref
  res <- list(
    call = call,
	pval = score,
	egm = egm
	)
  }else{
  res <- list(
    call = call,
	pval = score
	)
  }
  
  return(
   res
  )
}  


softmaxP <- function(p){
   expp <- exp(1-p)
   expp <- diag(1/rowSums(expp)) %*% expp
   expp
}