#' Canonical analysis of several sets with another set
#'
#' Relative proximities of several subsets of variables Yj with another set X. SUCCESSIVE SOLUTIONS
#'
#' The first solution calculates a standardized canonical component \code{cx[,1]} of x associated to ky
#' standardized components \code{cyi[,1]} of yi by maximizing \eqn{\sum_i \rho(cx[,1],cy_i[,1])^2}.
#' The second solution is obtained from the same criterion, with ky
#' orthogonality constraints for having \code{rho(cyi[,1],cyi[,2])=0} (that
#' implies \code{rho(cx[,1],cx[,2])=0)}.  For each of the 1+ky sets, the r
#' canonical components are 2 by 2 zero correlated.
#' The ky matrices (cx)'*cyi are triangular.
#' This function uses concor function.
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions
#'
#' @return A \code{list} with following components:
#' \item{cx}{a \code{n} times \code{r} matrix of the r canonical components of x}
#' \item{cy}{a \code{n.ky} times \code{r} matrix. The ky blocks cyi of the rows n*(i-1)+1 : n*i contain the r canonical components relative to Yi}
#' \item{rho2}{a \code{ky} times \code{r} matrix; each column k contains ky squared canonical correlations \eqn{\rho(cx[,k],cy_i[,k])^2}}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Hanafi & Lafosse (2001) Generalisation de la regression lineaire simple pour analyser la dependance de K ensembles de variables avec un K+1 eme.  Revue de Statistique Appliquee vol.49, n.1
#'
#' @examples
#'
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' ca <- concorcano(x,y,c(3,2,4),2)
#'
#' @export

concorcano  <-
  function(x, y, py, r) {
    # INITIALISATIONS

    n <- dim(x)[1]
    q <- dim(y)[2]
    if (sum(py) != q ) stop("py IS NOT SUITABLE")

    s <- svd(x);rk <- sum(s$d > max(dim(x))*s$d[1]*1e-8)
    P <- matrix(s$u[,1:rk]*sqrt(n),ncol=rk)

    ky <- length(py)
    ry <- matrix(0,1,ky)

    Py <- NULL
    cuy=c(0,cumsum(py))
    for (j in 1:ky) {
      s <- svd(y[,(cuy[j]+1):cuy[j+1]])
      ry[j] <- sum(s$d > max(c(n,py[j]))*s$d[1]*1e-8)
      Py <- cbind(Py,s$u[,1:ry[j]])
    }

    if (r > min(c(min(ry),rk,n))) stop("r IS TOO HIGH")
    Py <- matrix(Py,ncol=sum(ry))*sqrt(n)
    s <- concor(P,Py,ry,r);
    cumuly <- cumsum(c(0,ry))
    cy <- matrix(0,ky*n,r)
    for (k in 1:ky) {
      ay <- (cumuly[k]+1):cumuly[k+1]
      cy[(n*(k-1)+1):(n*k),] <- matrix(Py[,ay],ncol=ry[k])%*%s$v[ay,]
    }
    list(cx=P%*%s$u,cy=cy,rho2=s$cov2)
  }

