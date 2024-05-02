#' Relative links of several subsets of variables
#'
#' Relative links of several subsets of variables Yj with another set X. SUCCESSIVE SOLUTIONS
#'
#' The first solution calculates 1+kx normed vectors: the vector \code{u[:,1]} of Rp associated to the ky vectors \code{vi[:,1]}'s of Rqi,
#' by maximizing \eqn{\sum_i \mbox{cov}(x*u[,k],y_i*v_i[,k])^2}, with 1+ky norm constraints on the axes.
#' A component \code{(x)(u[,k])} is associated to ky partial components \code{(yi)(vi)[,k]} and to a global component \code{y*V[,k]}.
#' \eqn{\mbox{cov}((x)(u[,k]),(y)(V[,k]))^2 = \sum \mbox{cov}((x)(u[,k]),(y_i)(v_i[,k]))^2}.
#' \code{(y)(V[,k])} is a global component of the components \code{(yi)(vi[,k])}.
#' The second solution is obtained from the same criterion, but after replacing each yi by \eqn{y_i-(y_i)(v_i[,1])(v_i[,1]')}.
#' And so on for the successive solutions 1,2,...,r.  The biggest number of solutions may be r = inf(n, p, qi), when the (x')(yi')(s)
#' are supposed with full rank; then rmax = min(c(min(py),n,p)).  For a set of r solutions, the matrix u'X'YV is diagonal and the
#' matrices u'X'Yjvj are triangular (good partition of the link by the solutions).
#' concor.m is the svdcp.m function applied to the matrix x'y.
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions
#'
#' @return A \code{list} with following components:
#' \item{u}{A \code{p} times \code{r} matrix of axes in \code{Rp} relative to \code{x; (u^prime)(u) = Identity}}
#' \item{v}{A \code{q} times \code{r} matrix of \code{ky} row blocks \code{vi (qi x r)} of axes in \code{Rqi} relative to \code{yi; vi^prime*vi = Identity}}
#' \item{V}{A \code{q} times \code{r} matrix of axes in \code{Rq} relative to \code{y; Vprime*V = Identity}}
#' \item{cov2}{A \code{ky} times \code{r} matrix; each column \code{k} contains \code{ky} squared covariances \eqn{\mbox{cov}(x*u[,k],y_i*v_i[,k])^2}, the partial measures of link}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux: Definition de K+1 uples synthetiques. Revue de Statistique Appliquee vol.45,n.4.
#'
#' @examples
#'
#' # To make some "GPA" : so, by posing the compromise X = Y,
#' # "procrustes" rotations to the "compromise X" then are :
#' # Yj*(vj*u').
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' co <- concor(x,y,c(3,2,4),2)
#'
#' @export

concor  <-
  function(x, y, py, r) {
    # INITIALISATIONS
    n <- dim(x)[1]
    p <- dim(x)[2]
    q <- dim(y)[2]

    if (sum(py) != q ) stop("py IS NOT SUITABLE")
    if (r > min(c(min(py), q, n))) stop("r IS TOO HIGH")

    ky <- length(py)
    cri <- matrix(0, ky, r)
    cumuly = cumsum(c(0, py))
    u <- matrix(0, p, r)
    V <- matrix(0, q, r)
    v <- V

    for (i in 1:r) {
      s <- svd(t(x)%*%y)
      u[, i] <- s$u[, 1]
      V[, i] <- s$v[, 1]
      c1 = s$d[1]^2
      for (k in 1:ky) {
        ay <- (cumuly[k]+1):cumuly[k+1]
        ny <- t(V[ay, i])%*%V[ay, i]
        cri[k, i] <- ny*c1
        if (ny > 1e-8) {
          v[ay, i] <- V[ay, i]/sqrt(ny)
          y[, ay] <- y[, ay]-y[, ay]%*%(v[ay, i]%*%t(v[ay, i]))
        }
      }
    }
    list(u = u, v = v, V = V, cov2 = cri/n^2)
  }

