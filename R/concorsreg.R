#' Redundancy of sets yj by one set x
#'
#' Regression of several subsets of variables Yj by another set X. SUCCESSIVE SOLUTIONS
#'
#’ This function uses the concors function
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param px The row vector which contains the numbers pi, i = 1,...,kx, of the kx subsets xi of x : \eqn{\sum_i p_i}=sum(px)=p. px is the partition vector of x
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions
#'
#' @return A \code{list} with following components:
#' \item{cx}{a \code{n} times \code{r}matrix of the r explanatory components}
#' \item{v}{is a \eqn{q \times r}{q x r} matrix of ky row blocks \eqn{v_i} (\eqn{q_i \times r}{q_i x r}) of axes in Rqi relative to yi; \eqn{v_i'*v_i = \mbox{Id}}}
#' \item{varexp}{is a \eqn{ky \times r}{ky x r} matrix; each column k contains ky explained variances \eqn{\rho(cx[,k],y_i*v_i[,k])^2 \mbox{var}(y_i*v_i[,k])}}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#'
#' @examples
#'
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' crs <- concorsreg(x,c(2,3),y,c(3,2,4),2)
#' crs$varexp[1,1,]
#'
#' @export

concorsreg <-
  function(x,px,y,py,r) {
    if (sum(px) != dim(x)[2] | sum(py) != dim(y)[2] ) stop("px or py IS NOT SUITABLE")
    n <- dim(x)[1]
    kx <- length(px)
    rx <- matrix(0,1,kx)
    Px <- NULL
    cux = c(0,cumsum(px))
    for (j in 1:kx) {
      s <- svd(x[,(cux[j]+1):cux[j+1]])
      rx[j] <- sum(s$d > max(c(n,px[j]))*s$d[1]*1e-8)
      Px <- cbind(Px,s$u[,1:rx[j]]*sqrt(n))
    }
    if (r > min(c(min(py),min(rx),n))) stop("r IS TOO HIGH")

    cux <- c(0,cumsum(rx))
    Px <- matrix(Px,ncol = cux[kx+1])
    s <- concors(Px,rx,y,py,r)
    cx <- matrix(0,n*kx,r)
    for  (j in 1:kx) cx[((j-1)*n+1):(j*n),]<-matrix(Px[,(cux[j]+1):cux[j+1]],nrow = n)%*%s$u[(cux[j]+1):cux[j+1],]
    list(cx = cx,v = s$v,varexp = s$cov2)
  }

