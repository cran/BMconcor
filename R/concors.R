#' simultaneous  concorgm
#'
#' concorgm with the set of r solutions simultaneously optimized
#'
#' This function uses the svdbips function
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
#' \item{cov2}{a \code{ky} times \code{r} matrix; each column \code{k} contains \code{ky} squared covariances \eqn{\mbox{cov}(x*u[,k],y_i*v_i[,k])^2}, the partial measures of link}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux: Definition de K+1 uples synthetiques. Revue de Statistique Appliquee vol.45,n.4.
#'
#' @examples
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' cs <- concors(x,c(2,3),y,c(3,2,4),2)
#' cs$cov2[1,1,]
#'
#' @export

concors <-
  function(x,px,y,py,r) {
    if (sum(px) != dim(x)[2] | sum(py) != dim(y)[2] ) stop("px or py IS NOT SUITABLE")
    s <- svdbips(t(x)%*%y,px,py,r)
    list(u=s$u,v=s$v,cov2=s$s2/dim(x)[1]^2)
  }

