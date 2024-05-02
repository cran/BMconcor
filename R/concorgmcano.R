#' Canonical analysis of subsets Yj with subsets Xi
#'
#' Canonical analysis of subsets Yj with subsets Xi. Relative valuations by squared correlations of the proximities of subsets Xi with subsets Yj. SUCCESSIVE SOLUTIONS
#'
#' For the first solution, \eqn{sum_i sum_j \mbox{rho2}(cx_i[,1],cy_j[,1])} is the optimized
#' criterion. The other solutions are calculated from the same criterion, but with
#' orthogonalities for having two by two zero correlated the canonical components defined for
#' each xi, and also for those defined for each yj.  Each solution associates kx canonical
#' components to ky canonical components.  When kx =1 (px=p), take `concorcano` function
#' This function uses the concorgm function
#'
#' @param x are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column
#' @param y See \code{x}
#' @param px The row vector which contains the numbers pi, i=1,...,kx, of the kx subsets xi of x : \eqn{\sum_i p_i}=sum(px)=p. px is the partition vector of x
#' @param py The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.
#' @param r The number of wanted successive solutions rmax <= min(min(px),min(py),n)
#'
#' @return A \code{list} with following components:
#' \item{cx}{is a \code{n.kx} times  \code{r} matrix of kx row blocks cxi (n x r). Each row block contains r partial canonical components}
#' \item{cy}{is a  \code{n.ky} times  \code{r} matrix of ky row blocks cyj (n x r). Each row block contains r partial canonical components}
#' \item{rho2}{is a  \code{kx} time  \code{ky} tims  \code{r} array; for a fixed solution k, \code{rho2[,,k]} contains kxky squared correlations \eqn{rho2(cx[n*(i-1)+1:n*i,k],cy[n*(j-1)+1:n*j,k])}, simultaneously calculated between all the yj with all the xi}
#'
#' @author \enc{Lafosse, R.}{R. Lafosse}
#'
#' @references Kissita G., Analyse canonique generalisee avec tableau de reference generalisee. Thesis, Ceremade Paris 9 Dauphine (2003).
#'
#' @examples
#' x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
#' x <- scale(x);y <- scale(y)
#' cc <- concorgmcano(x,c(2,3),y,c(3,2,4),2)
#' cc$rho2[1,1,]
#'
#' @export

concorgmcano  <-
  function(x,px,y,py,r) {
    if (sum(px) != dim(x)[2] | sum(py) != dim(y)[2] ) stop("px or py IS NOT SUITABLE")
    n <- dim(x)[1]
    kx <- length(px)
    rx <- matrix(0,1,kx)
    Px <- NULL
    cux=c(0,cumsum(px))
    for (j in 1:kx) {
      s <- svd(x[,(cux[j]+1):cux[j+1]])
      rx[j] <- sum(s$d > max(c(n,px[j]))*s$d[1]*1e-8)
      Px <- cbind(Px,s$u[,1:rx[j]]*sqrt(n))
    }
    cux <- c(0,cumsum(rx))
    Px <- matrix(Px,nrow=n)
    ky <- length(py)
    ry <- matrix(0,1,ky)
    Py <- NULL
    cuy=c(0,cumsum(py))
    for (j in 1:ky) {
      s <- svd(y[,(cuy[j]+1):cuy[j+1]])
      ry[j] <- sum(s$d > max(c(n,py[j]))*s$d[1]*1e-8)
      Py <- cbind(Py,s$u[,1:ry[j]]*sqrt(n))
    }
    if (r > min(c(min(ry),min(rx),n))) stop("r IS TOO HIGH")
    cuy <- c(0,cumsum(ry))
    Py <- matrix(Py,nrow=n)
    s <- concorgm(Px,rx,Py,ry,r)
    cy <- matrix(0,n*ky,r)
    cx <- matrix(0,n*kx,r)

    for  (j in 1:kx) {
      cx[((j-1)*n+1):(j*n),] <- matrix(Px[,(cux[j]+1):cux[j+1]],nrow=n)%*%s$u[(cux[j]+1):cux[j+1],]
    }
    for  (j in 1:ky) cy[((j-1)*n+1):(j*n),] <- matrix(Py[,(cuy[j]+1):cuy[j+1]],nrow=n)%*%s$v[(cuy[j]+1):cuy[j+1],]

    list(cx=cx,cy=cy,rho2=s$cov2)
  }
