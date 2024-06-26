#' Redundancy of sets yj by one set x
#'
#' Regression of several subsets of variables Yj by another set X. SUCCESSIVE SOLUTIONS
#'
#’ The first solution calculates 1+ky normed vectors: the component
#’ cx[,1] in \eqn{R^n} associated to the ky vectors vi[,1]'s of \eqn{R^{q_i}}, by
#’ maximizing \eqn{varexp1=\sum_i \rho(cx[,1],y_i*v_i[,1])^2 \mbox{var}(y_i*v_i[,1]))}, with
#’ \eqn{1+ky} norm constraints. A explanatory component cx[,k] is associated to
#’ ky partial explained components yi*vi[,k] and also to a global explained
#’ component y*V[,k]. \eqn{\rho(cx[,k],y*V[,k])^2 \mbox{var}(y*V[,k])= \mbox{varexpk}}.  The
#’ total explained variance by the first solution is maximal.
#’ The second solution is obtained from the same criterion, but after
#’ replacing each yi by \eqn{y_i-y_i*v_i[,1]*v_i[,1]'}.  And so on for the
#’ successive solutions 1,2,...,r .  The biggest number of solutions may
#’ be \eqn{r=inf(n,p,q_i)}, when the matrices x'*yi are supposed with full rank.
#’ For a set of r solutions, the matrix (cx)'*y*V is diagonal : "on
#’ average", the explanatory component of one solution is only linked
#’ with the components explained by this explanatory, and is not linked
#’ with the explained components of the other solutions.  The matrices
#’ \eqn{(cx)'*y_j*v_j} are triangular : the explanatory component of one solution
#’ is not linked with each of the partial components explained in the
#’ following solutions.  The definition of the explanatory components
#’ depends on the partition vector py from the second solution.
#’ This function is using concor function  to the matrix x'y.
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions
#'
#' @return A \code{list} with following components:
#' \item{cx}{a \code{n} times \code{r}matrix of the r explanatory components}
#' \item{v}{is a \eqn{q \times r}{q x r} matrix of ky row blocks \eqn{v_i} (\eqn{q_i \times r}{q_i x r}) of axes in Rqi relative to yi; \eqn{v_i'*v_i = \mbox{Id}}}
#' \item{V}{is a \eqn{q \times r}{q x r} matrix of axes in Rq relative to y; \eqn{V'*V = \mbox{Id}}}
#' \item{varexp}{is a \eqn{ky \times r}{ky x r} matrix; each column k contains ky explained variances \eqn{\rho(cx[,k],y_i*v_i[,k])^2 \mbox{var}(y_i*v_i[,k])}}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux: Definition de K+1 uples synthetiques. Revue de Statistique Appliquee vol.45,n.4.
#' @references Chessel D. & Hanafi M. (1996) Analyses de la Co-inertie de K nuages de points.  Revue de Statistique Appliquee vol.44, n.2. (this ACOM analysis of one multiset is obtained by the command : concoreg(Y,Y,py,r))
#'
#' @examples
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' co <- concoreg(x,y,c(3,2,4),2)
#'
#' @export

concoreg  <-
  function(x, y, py, r) {
    n <- dim(x)[1]
    q <- dim(y)[2]
    if (sum(py) != q ) stop("py IS NOT SUITABLE")

    # rk here is a maximal value given for the rank of x.
    # you may modify the tolerance 1e-8
    s <- svd(x);rk <- sum(s$d > max(dim(x))*s$d[1]*1e-8)
    if (r > min(c(min(py),rk,n))) stop("r IS TOO HIGH")

    P=matrix(s$u[,1:rk]*sqrt(n),ncol=rk)
    s <- concor(P,y,py,r)
    list(cx=P%*%s$u,v=s$v,V=s$V,varexp=s$cov2)
  }

