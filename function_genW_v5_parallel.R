source("./function_OLS.R")

W_svd_v5_parallel <- function(df,SVDdecomp = FALSE){
	
	ngenes <- ncol(df$genes)
	ntranscripts <- ncol(df$transcripts)

 
  betayg_vec <- parApply(X = df$genes, MARGIN = 2, FUN = OLS_cross, y = df$y)
 
	
  betaxg_mx <- t(parApply(X = df$genes,MARGIN = 2,FUN = function(gene) {
      apply(X = df$transcripts, MARGIN=2, FUN=OLS_cross, x=gene)
    }))
  
  
	indexgrid <- as.matrix(expand.grid(
																			geneindex = seq(1,ngenes),
																			xindex = seq(1,ntranscripts)))
	indexgrid <- indexgrid[,c(2,1)]
	
	betay_xg_vec <- parApply(X = indexgrid,MARGIN = 1,FUN = function(cc){
									xg <- cbind(
															df$transcripts[,cc[1],drop=FALSE],
															df$genes[,cc[2],drop=FALSE])
									
									OLS_cross(xg,df$y)[1,]
									})
  
  betay_xg_mx <- matrix(betay_xg_vec,nrow = ngenes,ncol = ntranscripts,byrow = FALSE)
  
  
  
  betayg_mx <- matrix(rep(betayg_vec,each = ntranscripts),nrow = ngenes, ncol = ntranscripts,
											byrow = TRUE)
											
											
  w <- betayg_mx * betaxg_mx * betay_xg_mx
  
  w <- ifelse(w < 0, 0, ifelse(w > betayg_mx^2, 0,w))
  
  alphamx <- w / betayg_mx^2

  
  if(SVDdecomp == TRUE){
  	singulars <- svd(w)$d
  	u <- svd(w)$u
  	v <- svd(w)$v
  	return(list(w=w,alpha=alphamx,singulars=singulars,u=u,v=v))
  }else{
  	return(list(w=w,alpha=alphamx))
  }
}