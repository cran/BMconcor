#' SVD for a Column Partitioned matrix x
#'
#' SVD for a Column Partitioned matrix x. r global successive solutions
#'
#' The first solution calculates 1+kx normed vectors: the vector \code{u[,1]} of
#' \eqn{R^p} associated to the kx vectors \code{vi[,1]}'s of \eqn{R^{q_i}}. by maximizing
#' \eqn{\sum_i (u[,1]'*x_i*v_i[,1])^2}, with 1+kx norm constraints.  A
#' value \eqn{(u[,1]'*x_i*v_i[,1])^2} measures the relative link between
#' \eqn{R^p} and \eqn{R^{q_i}} associated to xi. It corresponds to a partial squared
#' singular value notion, since \eqn{\sum_i (u[,1]'*x_i*v_i[,1])^2=s^2},
#' where s is the usual first singular value of x.
#' The second solution is obtained from the same criterion, but after
#' replacing each xi by \code{xi-xi*vi[,1]*vi[,1]^prime}.  And so on for the
#' successive solutions 1,2,...,r .  The biggest number of solutions may
#' be r=inf(p,qi), when the xi's are supposed with full rank; then
#' \code{rmax=min([min(H),p])}.
#'
#' @param x a \code{p} times \code{q} matrix
#' @param H is a row vector which contains the numbers qh, h=1,...,ky, of the partition of x with ky column blocks : sum(qh)=q
#' @param r The number of wanted successive solutions
#'
#' @return A \code{list} with following components:
#' \item{u}{a \code{p} times \code{r} matrix of kx row blocks uk (pk x r); uk'*uk = Identity.}
#' \item{v}{a \code{q} times \code{r} matrix of \code{ky} row blocks \code{vi (qi x r)} of axes in \code{Rqi} relative to \code{yi; vi^prime*vi = Identity}}
#' \item{s}{a \code{kx} times \code{ky} times \code{r} array; with r fixed, each matrix contains kxky values \eqn{(u_h'*x_{kh}*v_k)^2}, the partial (squared) singular values relative to xkh.}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux: Definition de K+1 uples synthetiques. Revue de Statistique Appliquee vol.45,n.4.
#'
#' @examples
#' x <- matrix(runif(200),10,20)
#' s <- svdcp(x,c(5,5,10),1)
#' ss <- svd(x);ss$d[1]^2
#' sum(s$s2)
#'
#' @export

svdcp  <-
  function(x,H,r) {
    # Initialisations
    q <- dim(x)[2]
    if (sum(H) != q) stop("The value of `H` is not suitable")

    k <- length(H)
    s2 <- matrix(0,k,r)
    u <- matrix(0,dim(x)[1],r);
    v <- matrix(0,q,r);
    kx <- cumsum(c(0,H));

    # Calculus
    for (i in 1:r) {
      s <- svd(x)
      u[,i] <- s$u[,1]
      for (j in 1:k) {
        ax  <-  (kx[j]+1):kx[j+1]
        norm2  <-  t(s$v[ax,1]) %*% s$v[ax,1]
        s2[j,i] <- norm2 * s$d[1]^2
        if (s2[j,i] > 1e-8) {
          v[ax,i] <- s$v[ax,1]/sqrt(norm2)
          x[,ax]  <-  x[,ax]-x[,ax]%*%(v[ax,i]%*%t(v[ax,i]))
        }

      }
    }

    list(u=u,v=v,s2=s2)
  }

