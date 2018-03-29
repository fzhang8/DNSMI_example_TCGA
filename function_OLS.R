
OLS <- function(x,y){
	tx <- t(x)
	solve(tx %*% x) %*% tx %*% y
}

OLS_qr <- function(x,y){
	solve(qr(x), y)
}




OLS_cross <- function(x,y){

	out <- try(solve(crossprod(x,x),crossprod(x,y)))
	
	if(class(out) == "try-error"){
		return(matrix(0,ncol(x),1))
	}else{return(out)}
}