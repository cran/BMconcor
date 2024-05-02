#' Analyzing a set of partial links between Xi and Yj
#'
#' Analyzing a set of partial links between Xi and Yj, SUCCESSIVE SOLUTIONS
#'
#' The first solution calculates 1+kx normed vectors: the vector \code{u[:,1]} of Rp associated to the ky vectors \code{vi[:,1]}'s of Rqi,
#' by maximizing \code{sum(cov((x)(u[,k]),(y_i)(v_i[,k]))^2)}, with 1+ky norm constraints on the axes.
#' A component \code{(x)(u[,k])} is associated to ky partial components \code{(yi)(vi)[,k]} and to a global component \code{y*V[,k]}.
#' \code{cov((x)(u[,k]),(y)(V[,k]))^2 = sum(cov((x)(u[,k]),(y_i)(v_i[,k]))^2)(y)(V[,k])} is a global component of the components \code{(yi)(vi[,k])}.
#' The second solution is obtained from the same criterion, but after replacing each yi by \eqn{y_i-(y_i)(v_i[,1])(v_i[,1]')}.
#' And so on for the successive solutions 1,2,...,r.  The biggest number of solutions may be r=inf(n, p, qi), when the (x')(yi')(s)
#' are supposed with full rank; then rmax=min(c(min(py),n,p)).  For a set of r solutions, the matrix u'X'YV is diagonal and the
#' matrices u'X'Yjvj are triangular (good partition of the link by the solutions).
#' concor.m is the svdcp.m function applied to the matrix x'y.
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param px A row vector which contains the numbers pi, i=1,...,kx, of the kx subsets xi of x : sum(pi)=sum(px)=p. px is the partition vector of x
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions rmax <= min(min(px),min(py),n)
#'
#' @return A \code{list} with following components:
#' \item{u}{a \code{p} times \code{r} matrix of axes in \code{Rp} relative to \code{x; u^prime*u = Identity}}
#' \item{v}{a \code{q} times \code{r} matrix of \code{ky} row blocks \code{vi (qi x r)} of axes in \code{Rqi} relative to \code{yi; vi^prime*vi = Identity}}
#' \item{cov2}{a \code{ky} times \code{r} matrix; each column \code{k} contains \code{ky} squared covariances \eqn{\mbox{cov}((x)(u[,k]),(y_i)(v_i[,k]))^2}, the partial measures of link}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Kissita, Cazes, Hanafi & Lafosse (2004) Deux methodes d'analyse factorielle du lien entre deux tableaux de variables partitionn?es. Revue de Statistique Appliqu?e, Vol 52, n. 3, 73-92.
#'
#' @examples
#'
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' cg <- concorgm(x,c(2,3),y,c(3,2,4),2)
#' cg$cov2[1,1,]
#'
#' @export

concorgm  <-
  function(x,px,y,py,r) {
    if (sum(px) != dim(x)[2] | sum(py) != dim(y)[2] ) stop("px or py IS NOT SUITABLE")
    s <- svdbip(t(x)%*%y,px,py,r)
    list(u=s$u,v=s$v,cov2=s$s2/dim(x)[1]^2)
  }

