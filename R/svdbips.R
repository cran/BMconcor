#' SVD for bipartitioned matrix x
#'
#' SVD for bipartitioned matrix x. SIMULTANEOUS SOLUTIONS. ("simultaneous svdbip")
#'
#' One set of r solutions is calculated by maximizing \eqn{\sum_i \sum_k \sum_h
#' (u_k[,i]'*x_{kh}*v_h[,i])^2}, with kx+ky orthonormality constraints (for
#' each uk and each vh).  For each fixed r value, the solution is totally
#' new (does'nt consist to complete a previous calculus of one set of r-1
#' solutions).  \code{rmax=min([min(K),min(H)])}.  When r=1, it is svdbip (thus
#' it is svdcp when r=1 and kx=1).
#' Convergence of algorithm may be not global. So the below proposed
#' initialisation of the algorithm may be not very suitable for some data
#' sets.  Several different random initialisations with normed vectors
#' might be considered and the best result then choosen....
#'
#' @param x a \code{p} times \code{q} matrix
#' @param K is a row vector which contains the numbers pk, k=1,...,kx, of the partition of x with kx row blocks : \code{sum(pk)=p}
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
#' @references Lafosse R. & Ten Berge J. A simultaneous CONCOR method for the analysis of two partitioned matrices. submitted.
#'
#' @examples
#' x <- matrix(runif(200),10,20)
#' s1 <- svdbip(x,c(3,4,3),c(5,5,10),2);sum(sum(sum(s1$s2)))
#' ss <- svdbips(x,c(3,4,3),c(5,5,10),2);sum(sum(sum(ss$s2)))
#'
#' @export

svdbips <-
  function(x,K,H,r) {
    # INITIALISATIONS
    p <- dim(x)[1]
    q <- dim(x)[2]
    if (sum(H) != q | sum(K) != p) stop("K or H IS NOT SUITABLE")
    if (r > min(c(K,H))) stop("r IS NOT SUITABLE")
    M <- length(K)
    N <- length(H)
    u <- matrix(0,p,r)
    v <- matrix(0,q,r)
    ck <- cumsum(c(0,K))
    ch <- cumsum(c(0,H))
    s2 <- array(0,c(M,N,r))

    #PROPOSED INITIALISATION OF THE ALGORITHM with u and v
    for (i in 1:M) {
      ak <- (ck[i]+1):ck[i+1]
      s <- svd(matrix(x[ak,],nrow=length(ak)))
      u[ak,]<-s$u[,1:r]
    }

    for (j in 1:N) {
      ah <- (ch[j]+1):ch[j+1]
      s <- svd(matrix(x[,ah],ncol=length(ah)))
      v[ah,]<-s$v[,1:r]
    }
    cc <- 2;cc1 <- 0

    #ALGORITHM
    while (abs(cc-cc1) > 1e-8) {
      #aa and bb are converging to the optimized criterion
      aa=0;bb=0;
      cc1=cc;
      A <- matrix(0,p,r)
      B <- matrix(0,r,q)

      for (i in 1:M) {
        ak <- (ck[i]+1):ck[i+1]
        for (j in 1:N) {
          ah <- (ch[j]+1):ch[j+1]
          d <- diag(t(u[ak,])%*%x[ak,ah]%*%v[ah,]);l <- length(d)
          A[ak,]<-A[ak,]+matrix(x[ak,ah],nrow=K[i])%*%v[ah,]%*%diag(d,nrow=l)
        }
        s <- svd(A[ak,]);u[ak,]<-s$u[,1:r]%*%t(s$v)
        aa <- aa+sum(s$d)
      }

      for (j in 1:N) {
        ah <- (ch[j]+1):ch[j+1]
        for (i in 1:M) {
          ak <- (ck[i]+1):ck[i+1]
          d <- diag(t(u[ak,])%*%x[ak,ah]%*%v[ah,]);l <- length(d)
          B[,ah]<-B[,ah]+diag(d,nrow=l)%*%t(u[ak,])%*%x[ak,ah]
        }
        s <- svd(t(B[,ah]));v[ah,]<-s$u[,1:r]%*%t(s$v)
        bb <- bb+sum(s$d)
      }

      cc <- (sqrt(aa)+sqrt(bb))/2
    }

    for (k in 1:r) {
      for (i in 1:M) {
        ak <- (ck[i]+1):ck[i+1]
        for (j in 1:N) {
          ah <- (ch[j]+1):ch[j+1]
          s2[i,j,k]<-(t(u[ak,k])%*%x[ak,ah]%*%v[ah,k])^2
        }
      }
    }
    list(u=u,v=v,s2=s2)
  }

