#' SVD for one bipartitioned matrix x
#'
#' SVD for bipartitioned matrix x. r successive Solutions
#'
#' The first solution calculates kx+ky normed vectors: kx vectors \code{uk[:,1]}
#' of \eqn{R^{p_k}} associated to ky vectors \code{vh[:,1]}'s of \eqn{R^{q_h}}, by
#' maximizing \eqn{\sum_k \sum_h (u_k[:,1]^prime*x_{kh}*v_h[:,1])^2}, with kx+ky
#' norm constraints. A value \eqn{(u_k[,1]^prime*x_{kh}*v_h[,1])^2} measures the
#' relative link between \eqn{R^{p_k}} and \eqn{R^{q_h}} associated to the block xkh.
#' The second solution is obtained from the same criterion, but after
#' replacing each xhk by xkh-xkh*vh*vh'-uk*uk'xkh+uk*uk'xkh*vh*vh'.  And
#' so on for the successive solutions 1,2,...,r .  The biggest number of
#' solutions may be r=inf(pk,qh), when the xkh's are supposed with full
#' rank; then \code{rmax=min([min(K),min(H)])}.
#' When K=p (or H=q, with t(x)), svdcp function is better.  When H=q and
#' K=p, it is the usual svd (with squared singular values).
#' Convergence of algorithm may be not global. So the below proposed
#' initialisation of the algorithm may be not very suitable for some data
#' sets.  Several different random initialisations with normed vectors
#' might be considered and the best result then choosen.
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
#' @references Kissita G., Cazes P., Hanafi M. & Lafosse (2004) Deux methodes d'analyse factorielle du lien entre deux tableaux de variables partitiones. Revue de Statistique Appliquee.
#'
#' @examples
#'
#' x <- matrix(runif(200),10,20)
#' s <- svdbip(x,c(3,4,3),c(5,15),3)
#'
#' @export

svdbip <-
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
    A <- matrix(0,p,N)
    B <- matrix(0,q,M)
    s2 <- array(0,c(M,N,r))

    # ALGORITHM
    for (k in 1:r) {

      #PROPOSED INITIALISATION OF THE ALGORITHM with u
      for (i in 1:M) {
        ak <- (ck[i]+1):ck[i+1]
        s <- svd(matrix(x[ak,],nrow=length(ak)))
        u[ak,k]<-s$u[,1]
      }


      cc <- s$d[1];cc1 <- 0;
      #comp <- 0;

      while (abs(cc-cc1) > 1e-8) {
        #aa^2 and bb^2 are converging to the optimized criterion
        aa <- 0;bb <- 0;
        cc1 <- cc;

        #comp <- comp+1;

        for (j in 1:N) {
          ah <- (ch[j]+1):ch[j+1]
          for (i in 1:M) {
            ak <- (ck[i]+1):ck[i+1]
            B[ah,i]<-t(matrix(x[ak,ah],nrow=length(ak)))%*%u[ak,k]
          }
          s <- svd(matrix(B[ah,],nrow=length(ah)))
          if (s$d[1] > 1e-8) { v[ah,k]<-s$u[,1]; aa <- aa+s$d[1] }
        }


        for (i in 1:M) {
          ak <- (ck[i]+1):ck[i+1]
          for (j in 1:N) {
            ah <- (ch[j]+1):ch[j+1]
            A[ak,j]<- matrix(x[ak,ah],nrow=length(ak))%*%v[ah,k]
          }
          s <- svd(matrix(A[ak,],nrow=length(ak)))
          if (s$d[1] > 1e-8) {u[ak,k]<-s$u[,1];bb <- bb+s$d[1]}
        }

        cc <- (aa+bb)/2
      }

      for (i in 1:M) {
        ak <- (ck[i]+1):ck[i+1]
        for (j in 1:N) {
          ah <- (ch[j]+1):ch[j+1]
          c <- t(u[ak,k])%*%x[ak,ah]%*%v[ah,k]
          x[ak,ah] <- x[ak,ah]-u[ak,k]%*%t(u[ak,k])%*%x[ak,ah]-x[ak,ah]%*%(v[ah,k]%*%t(v[ah,k]))+u[ak,k]%*%c%*%t(v[ah,k])
          s2[i,j,k] <- c^2
        }
      }
      #comp
    }
    list(u=u,v=v,s2=s2)
  }

