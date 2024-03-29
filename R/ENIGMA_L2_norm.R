#' @title ENIGMA maximum L2 norm version
#'
#' @description maximum L2 norm version is the default version of ENIGMA, which is the regularized weighted matrix completion to constraint the maximum L2 norm of inferred cell type-specific gene expression matrix.
#' 
#' @param object ENIGMA object
#' @param alpha
#' ENIGMA is a multi-objective optimization problem involve two object function, the distance function between observed bulk RNA-seq and reconstitute RNA-seq generated by weighted combination of CSE, and the distance function beween average CSE expression and cell type reference matrix. The alpha is used to determine weights of these two objects. If the alpha gets larger, the optimization attach greater importance on the the first object. Default: 0.5
#'
#' @param tao_k
#' The step size of each round of gradient decent. Default: 0.01
#'
#' @param beta
#' The regularization parameter to penalize the weight (deconvoluted expression) matrices from being too large. Default: 0.1
#'
#' @param epsilon
#' Determine the stop condition in CSE updating. Default: 0.001
#'
#' @param max.iter
#' The maximum number of iterations. Default: 1000
#'
#' @param verbose
#' Rreturn the information after each step of processing. Default: TRUE
#'
#' @param pos
#' Set all entries in CSE is positive. Default: TRUE
#'
#' @param calibrate
#' calibrate the inferred CSE into input bulk gene expression scale. Default: TRUE
#'
#' @param Norm.method
#' Method used to perform normalization. User could choose PC, frac or quantile
#' 
#' @param preprocess
#' Method used to perform variance stablization preprocessing. User could choose none, sqrt or log
#'
#' @param loss_his
#' save the loss value of each round of iteration.
#'
#' @param model_tracker
#' save the model in returned object
#' 
#' @param model_name
#' name of the model
#'
#' @param X_int
#' initialization for CSE profiles, an array object with three dimensions (the number of genes * the number of samples * the number of cell types), if user input a matrix (the number of genes * the number of samples), each cell type would be assigned the same start matrix.
#'
#' @param force_normalize
#' when alpha >= 0.9 or profile matrix is not generated from S-mode batch effect correction, ENIGMA would not perform normalization, if user still want to perform normalization, set force_normalize=TRUE. Default: FALSE
#'
#' @return ENIGMA object where object@result_CSE contains the inferred CSE profile, object@result_CSE_normalized would contains normalized CSE profile, object@loss_his would contains the loss values of object functions during model training. If model_tracker = TRUE, then above results would be saved in the object@model.
#'
#'
#' @examples
#' \dontrun{
#' egm = ENIGMA_l2max_norm(egm,model_tracker = TRUE, Norm.method="quantile")
#' egm@result_CSE
#' egm@result_CSE_normalized
#' }
#'
#'
#' @export
ENIGMA_L2_max_norm <- function(object, alpha=0.5, tao_k=0.01, beta=0.1, epsilon=0.001, max.iter=1000,verbose=FALSE, pos=TRUE,calibrate=TRUE, Norm.method = "frac",preprocess = "sqrt",loss_his=TRUE,model_tracker=FALSE,model_name = NULL,X_int=NULL, Normalize =  TRUE){
    suppressPackageStartupMessages(require("scater"))
	suppressPackageStartupMessages(require("preprocessCore"))
	
	###Create a model assay
	if(model_tracker){
	if(is.null(model_name)){
	  model_name = paste("maximum_L2_norm_model_",date(),"_trained",sep="")
	}
	  basic_infor = data.frame(alpha = alpha,beta = beta, step_size = tao_k, epsilon = epsilon,max_iter = max.iter,calibrate = calibrate,Normalize_method = Norm.method,preprocess = preprocess,pos=pos)
	  object@model[[model_name]] <- list(basic_infor = basic_infor)
	}
	
	if ( !(preprocess %in% c("none", "sqrt","log")) | (length(preprocess) != 1) ) {
        stop("Invalid data transformation method type. Please input 'none','sqrt' or 'log'. ")
    }
	if ( !(Norm.method %in% c("PC", "frac","quantile")) | (length(Norm.method) != 1) ) {
        stop("Invalid normalization method type. Please input 'PC','frac' or 'quantile'. ")
    }
	
    X = object@bulk
    theta = object@result_cell_proportion
	R = object@ref
	
    # unify geneid between X and R
    geneid = intersect( rownames(X), rownames(R) )
    X = X[geneid,]
    R = R[geneid,]
  
    ## renormalization
	geneID <- rownames(X)
	sampleID <- colnames(X)
	ctID <- colnames(R)
	X <- X %*% diag(10^5/colSums(X))
	R <- R %*% diag(10^5/colSums(R))
	rownames(X) <- rownames(R) <- geneID
	colnames(X) <- sampleID
	colnames(R) <- ctID
    
    if(pre.process == "log"){
	 X <- log2(X+1)
	 R <- log2(R+1)
	}
	if(pre.process == "sqrt"){
	 X <- sqrt(X)
	 R <- sqrt(R)
	}

    # initialize the CSE
    P_old = array(0,
                  dim = c( nrow(X),
                           ncol(X),
                           ncol(theta)),
                  dimnames = list( rownames(X),
                                   colnames(X),
                                   colnames(theta))
    )
    X_int_m = array(0,
              dim = c( nrow(X),
                       ncol(X),
                       ncol(theta)),
              dimnames = list( rownames(X),
                               colnames(X),
                               colnames(theta))
    )
    if(is.null(X_int) == FALSE){
       if(length(dim(X_int)) == 2){
           for(i in 1:ncol(theta)){
           X_int_m[,,i] = X_int
           }
        }
       if(length(dim(X_int)) == 3){
          for(i in 1:ncol(theta)){
           X_int_m[,,i] = X_int[,,i]
          }
        }
    }
    for(i in 1:ncol(theta)){
        if(is.null(X_int)){P_old[,,i] <- X}else{P_old[,,i] <- X_int_m[,,i]}
    }
    rm(X_int, X_int_m);gc()
    ###update iteractively
    P_old_new <- P_old
    mask_entry <- matrix(1,nrow = nrow(X), ncol = ncol(X)); mask_entry[X==0] <- 0
    cat(date(), 'Optimizing cell type specific expression profile... \n')
            iter.exp <- 0
			loss <- NULL
			DisList <- NULL
            repeat{
                ratio <- NULL
                dP <- derive_P2(X, theta,P_old,R,alpha,mask_entry)
                for(i in 1:ncol(theta)){
                    P_hat <- proximalpoint(P_old[,,i], tao_k,dP[,,i],beta*10^5)
                    P_old_new[,,i] <- P_hat
                    ratio <- c(ratio, sum( (P_hat-P_old[,,i])^2 ))
                }
                if(verbose) writeLines( sprintf("   L2 distance ranges from: %f - %f", min(ratio), max(ratio) ) )
				r <- sub_loss(X, P_old, theta, alpha,beta,R)
                if(loss_his) loss<- rbind(loss,c(r$part1,r$part2,r$part3))
				DisList <- c(DisList,max(ratio))
				if(iter.exp>10){
				converge_score <- DisList[(length(DisList)-9):length(DisList)] - DisList[(length(DisList)-10):(length(DisList)-1)]
				if(sum(converge_score)>2){
				  stop(paste('Convergence Error! please set a higher value for beta or set a smaller value for learning rate (tao_k)\n'
           ,'Suggest set beta = ',beta*2,' or tao_k = ',tao_k/2,sep=""))
				}
				}
                if(max(ratio) < epsilon||iter.exp >= max.iter){break}else{

	    #update P_old
                    P_old <- P_old_new
                    iter.exp <- iter.exp + 1
                }
            }
			
        ## Return the Loss
        loss_new.obj <- sub_loss(X, P_old, theta, alpha,beta,R)
        if(verbose) writeLines( paste("Total loss: ", loss_new.obj$val ,sep="") )
        if(verbose) writeLines( paste("part1:",loss_new.obj$part1," part2:",loss_new.obj$part2," part3:",loss_new.obj$part3,sep="") )
        ### optimized theta
        ### optimize theta
        ### take the gradient of all theta and running gradient decent
        if(pos){P_old[P_old<0] <- 0}
    X_k_norm <- X_k <- P_old
	if(pos&verbose&calibrate) writeLines("calibration...")
	for(k in 1:dim(X_k)[3]){
	if(preprocess == "sqrt") X_k[,,k] <- X_k[,,k]^2
	if(preprocess == "log") X_k[,,k] <- 2^X_k[,,k] - 1
	if(pos&calibrate){
	X_k[,,k] <- X_k[,,k] * (mean(colSums(object@bulk))/mean(colSums(X_k[,,k])))
	}
	}
	
	if(Normalize){
	if(verbose) cat("Perform Normalization...")
	#X_k_norm <- X_k <- P_old
	if(Norm.method == "PC"){
	for(k in 1:dim(X_k)[3]){
	   exp <- X_k[,,k]
	   scale_x <- function(x){
	   if(var(x)==0){x <- x - mean(x)}else{x <- scale(x)}
	   x
	   }
	exp.scale <- t(apply(exp,1,scale_x))
	###chose the PC with the highest correlation with cell type fractions
	d <- sqrt(svd(exp.scale)$d)
	d <- d / sum(d)
	prob_d <- NULL;for(i in 1:length(d)) prob_d <- c(prob_d, sum(d[1:i]))
	PC <- svd(exp.scale)$v[,1:which(prob_d>0.8)[1]]
	pc_cor <- apply(PC,2,function(x){cor(x,theta[,k],method="sp")})
	if(max(abs(pc_cor))<0.9){
    warning("cannot found cell type proportion related principal component, using frac-based normalization automatically")
	Norm.method = "frac"
	}else{
	PC <- PC[,which.max(abs(pc_cor))]
	the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
	exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
	X_k_norm[,,k] <- exp.norm
	}
	}
	}
        
	if(Norm.method == "frac"){
	for(k in 1:dim(X_k)[3]){
	   exp <- X_k[,,k]
	   PC <- theta[,k]
       the <- (exp %*% as.matrix(PC) - length(PC) * mean(PC) * rowMeans(exp)) / (sum(PC^2) - length(PC)*mean(PC)^2)
       exp.norm <- exp - as.matrix(the) %*% t(as.matrix(PC))
	   X_k_norm[,,k] <- exp.norm
	}
	}
	
	if(Norm.method == "quantile"){
	for(k in 1:dim(X_k)[3]){
	   exp <- X_k[,,k]
	   exp.norm <- normalize.quantiles(exp)
	   rownames(exp.norm) <- rownames(exp)
	   colnames(exp.norm) <- colnames(exp)
	   X_k_norm[,,k] <- exp.norm
	}
	}
	
	object@result_CSE_normalized = res2sce(X_k_norm)
	if(model_tracker){
	  object@model[[model_name]]$result_CSE_normalized = res2sce(X_k_norm)
	}
	
	}else{
	object@result_CSE_normalized = res2sce(X_k)
	if(model_tracker){
	  object@model[[model_name]]$result_CSE_normalized = object@result_CSE_normalized
	}
	}
	
	writeLines( paste("Converge in ",iter.exp," steps",sep="") )
	# return cell type specific gene expression profile
    object@result_CSE = res2sce(X_k)
	
	##loading loss history
	if(loss_his) object@loss_his = loss
	if(model_tracker){
	   if(loss_his){
	   object@model[[model_name]]$loss_his = loss
	   }else{
	   object@model[[model_name]]$loss_his = NULL
	   }
	   object@model[[model_name]]$result_CSE = res2sce(X_k)
	   
	   ### import the model name and model type
	   if(nrow(object@model_name)==0){
	   m = t(as.matrix(c(model_name, "maximum L2 norm model")))
	   m = as.data.frame(m)
	   colnames(m) = c("Model Name","Model Type")
	   rownames(m) = paste0("model",1:nrow(m))
	   object@model_name = m
	   }else{
	   object@model_name = subset(object@model_name,`Model Name` %in% model_name == FALSE)
	   
	   if(nrow(object@model_name)==0){
	   m = t(as.matrix(c(model_name, "maximum L2 norm model")))
	   m = as.data.frame(m)
	   colnames(m) = c("Model Name","Model Type")
	   rownames(m) = paste0("model",1:nrow(m))
	   object@model_name = m
	   }else{
	   m = rbind(as.matrix(object@model_name),t(as.matrix(c(model_name, "maximum L2 norm model"))))
	   m = as.data.frame(m)
	   colnames(m) = c("Model Name","Model Type")
	   rownames(m) = paste0("model",1:nrow(m))
	   object@model_name = m
	   }
	   }
	}
	if(verbose) cat(date(),'Done... \n')
    return(object)
}

sub_loss <- function(X, P_old, theta, alpha,beta,R){
    # X: Bulk gene expression dataset (g*n)
    # P_old: cell type specific gene expression profile (g*n*p)
    # theta: cell type ratio for each samples (n*p)
    # alpha: constraint parameters of the similarity between each estimated cell type specific expression and reference profile, constant
    # beta:  constraint parameters of the smoothness of gene expression, constant
    # R: reference profile (g*p)

    part1 <- 0
    for(i in 1:ncol(theta)){
        part1 <- part1+P_old[,,i]%*%diag(theta[,i])
    }
    part1 <- part1
    # part1 <- norm((X-part1),"F")^2
    part1 <- sum( (X-part1)^2 )

    part2 <- 0
    for(i in 1:ncol(R)){
        # part2 <- part2 + alpha*norm((P_old[,,i]-ref),"F")^2
        part2 <- part2 + sum( (rowMeans(P_old[,,i])-R[,i])^2 )
    }

    part3 <- 0
    for(i in 1:ncol(R)){
        # norm <- apply(P_old[,,i],2,norm,"2")
        part3 <- part3 + max( colSums(P_old[,,i]^2) )
    }

    res <- list()
    val <- part1+part2+beta*part3
    res$val <- val
    res$part1 <- part1*(alpha/2)
    res$part2 <- part2*((1-alpha)/2)
    res$part3 <- part3*(beta/2)
    res
}

squash <- function(V, beta){
    ## squash: calculate the optimal solution of the formula: X=argmin{ (||X-V||_F)^2 + beta*||X||_2_max }
    n <- NULL
    for(i in 1:nrow(V)){
        n <- c(n, sqrt( sum(V[i,]^2) ) )
    }
    pi <- order(n,decreasing=TRUE)
    s <- NULL
    for(i in 1:length(pi)){
        s <- c(s, sum(n[pi[1:i]]))
    }
    q <- max(which(n[pi]>=s/(c(1:length(s))+beta)))
    tao <- s[q]/(q+beta)

    for(i in 1:q){
        V[pi[i],] <- tao*V[pi[i],]/sqrt( sum(V[pi[i],]^2) )
    }

    V
}

proximalpoint <- function(P, tao_k,dP,beta){
    # X: Bulk gene expression dataset (g*n)
    # P_old: cell type specific gene expression profile (g*n*p)
    # theta: cell type ratio for each samples (n*p)
    # alpha: constraint parameters of the similarity between each estimated cell type specific expression and reference profile, constant
    # beta:  constraint parameters of the smoothness of gene expression, constant
    # R: reference profile (g*p)
    # P: the ith cell type specific gene expression profile needs to be undated
    # tao_k: gradient size
    # dP: gradient of matrix P
    # scale_alpha: the parameters for inequality decision
    # beta:  constraint parameters of the smoothness of gene expression, constant
    # cell_type_index: optimize which type of cells
    # gamma: the parameters for inequality decision

    P_hat <- t(squash(t(P-tao_k*dP),tao_k*beta))
    ##update P matrix
    return(P_hat)
}



derive_P2 <- function(X, theta, P_old,R,alpha,mask_entry){
  ## P_old: a tensor variable with three dimensions
  ## theta: the cell type proportions variable
  ## cell_type_index: optimize which type of cells
  ## R: reference matrix
  dP1 <- dP2 <- array(0,
                      dim = c( nrow(X),
                               ncol(X),
                               ncol(theta)),
                      dimnames = list( rownames(X),
                                       colnames(X),
                                       colnames(theta))
  )
  for(cell_type_index in 1:ncol(theta)){
    R.m <- as.matrix(R[,cell_type_index])

    cell_type_seq <- c(1:ncol(theta))
    cell_type_seq <- cell_type_seq[cell_type_seq!=cell_type_index]

    X_summary = Reduce("+",
                       lapply(cell_type_seq, function(i) P_old[,,i]%*%diag(theta[,i]) )
    )
    X_summary <- X-X_summary

    dP1[,,cell_type_index] <- 2*(P_old[,,cell_type_index]%*%diag(theta[,cell_type_index]) - X_summary)%*%diag(theta[,cell_type_index])
    dP2[,,cell_type_index] <- 2*(as.matrix(rowMeans(P_old[,,cell_type_index]))-R.m)%*%t(as.matrix(rep((1/ncol(dP2[,,cell_type_index])),ncol(dP2[,,cell_type_index]))))
	dP1[,,cell_type_index] <- dP1[,,cell_type_index]*mask_entry
	dP2[,,cell_type_index] <- dP2[,,cell_type_index]*mask_entry
  }
  dP1 = dP1 / sqrt( sum( dP1^2 ) ) * 1e5
  dP2 = dP2 / sqrt( sum( dP2^2 ) ) * 1e5

  w1 <- alpha
  w2 <- 1-w1

  dP <- dP1*as.numeric(w1) + dP2*as.numeric(w2)
  return(dP)
}