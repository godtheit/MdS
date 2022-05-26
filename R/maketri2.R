
maketri2 <- function(y,D,R,k) {
  m = nrow(D)
  n = ncol(D)

  a = .C(C_maketri2,
         y=as.double(y),
         D=as.double(D),
         R=as.double(R),
         m=as.integer(m),
         n=as.integer(n),
         k=as.integer(k),
         PACKAGE="MdS2")

  y = a$y
  D = matrix(a$D,nrow=m)
  R = matrix(a$R,nrow=m)

  return(list(y=y,D=D,R=R))
}
