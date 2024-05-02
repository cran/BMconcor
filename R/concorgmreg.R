#' Regression of subsets Yj by subsets Xi
#'
#' Regression of subsets Yj by subsets Xi for comparing all the explanatory-explained pairs (Xi,Yj). SUCCESSIVE SOLUTIONS
#'
#' For the first solution, \eqn{\sum_i \sum_j \mbox{rho2}(cx_i[,1],y_j*v_j[,1])
#' \mbox{var}(y_j*v_j[,1])} is the optimized criterion. The second solution is
#' calculated from the same criterion, but with \eqn{y_j-y_j*v_j[,1]*v_j[,1]'}
#' instead of the matrices yj and with orthogonalities for having two by
#' two zero correlated the explanatory components defined for each matrix
#' xi. And so on for the other solutions. One solution k associates kx
#' explanatory components (in \code{cx[,k]}) to ky explained components. When
#' kx =1 (px = p), take concoreg function
#' This function uses the concorgm function
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param px A row vector which contains the numbers pi, i = 1,...,kx, of the kx subsets xi of x : sum(pi)=sum(px)=p. px is the partition vector of x
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions
#'
#' @return A \code{list} with following components:
#' \item{cx}{a \code{n} times \code{r}matrix of the r explanatory components}
#' \item{v}{is a \eqn{q \times r}{q x r} matrix of ky row blocks \eqn{v_i} (\eqn{q_i \times r}{q_i x r}) of axes in Rqi relative to yi; \eqn{v_i'*v_i = \mbox{Id}}}
#' \item{varexp}{is a kx x ky x r array; for a fixed solution k, the matrix \code{varexp[,,k]} contains kxky explained variances obtained by a simultaneous regression of all the yj by all the xi, so the values \eqn{\mbox{rho2}(cx[n*(i-1)+1:n*i,k],y_j*v_j[,k]) var(y_j*v_j[,k])}}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Hanafi & Lafosse (2004) Regression of a multi-set by another based on an extension of the SVD. COMPSTAT'2004 Symposium
#'
#' @examples
#'
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' cr <- concorgmreg(x,c(2,3),y,c(3,2,4),2)
#' cr$varexp[1,1,]
#'
#' @export

concorgmreg  <-
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
      Px <- cbind(Px,s$u[,1:rx[j]])
    }
    Px <- Px*sqrt(n)
    if (r > min(c(min(py),min(rx),n))) stop("r IS TOO HIGH")

    cux <- c(0,cumsum(rx))
    Px <- matrix(Px,ncol = cux[kx+1])
    s <- concorgm(Px,rx,y,py,r)
    cx <- matrix(0,n*kx,r)
    for  (j in 1:kx) cx[((j-1)*n+1):(j*n),] <- matrix(Px[,(cux[j]+1):cux[j+1]],nrow = n)%*%s$u[(cux[j]+1):cux[j+1],]
    list(cx = cx,v = s$v,varexp = s$cov2)
  }

