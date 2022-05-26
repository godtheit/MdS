Backsolve <- function(A,b,q) {
  m = nrow(A)
  n = ncol(A)
  z = backsolve(A[Seq(1,n-q),Seq(q+1,n),drop=FALSE],
                b[Seq(q+1,n)],transpose=TRUE)
  return(c(z,rep(0,m-n+q)))
}
