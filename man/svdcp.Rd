% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svdcp.R
\name{svdcp}
\alias{svdcp}
\title{SVD for a Column Partitioned matrix x}
\usage{
svdcp(x, H, r)
}
\arguments{
\item{x}{a \code{p} times \code{q} matrix}

\item{H}{is a row vector which contains the numbers qh, h=1,...,ky, of the partition of x with ky column blocks : sum(qh)=q}

\item{r}{The number of wanted successive solutions}
}
\value{
A \code{list} with following components:
\item{u}{a \code{p} times \code{r} matrix of kx row blocks uk (pk x r); uk'*uk = Identity.}
\item{v}{a \code{q} times \code{r} matrix of \code{ky} row blocks \code{vi (qi x r)} of axes in \code{Rqi} relative to \code{yi; vi^prime*vi = Identity}}
\item{s}{a \code{kx} times \code{ky} times \code{r} array; with r fixed, each matrix contains kxky values \eqn{(u_h'*x_{kh}*v_k)^2}, the partial (squared) singular values relative to xkh.}
}
\description{
SVD for a Column Partitioned matrix x. r global successive solutions
}
\details{
The first solution calculates 1+kx normed vectors: the vector \code{u[,1]} of
\eqn{R^p} associated to the kx vectors \code{vi[,1]}'s of \eqn{R^{q_i}}. by maximizing
\eqn{\sum_i (u[,1]'*x_i*v_i[,1])^2}, with 1+kx norm constraints.  A
value \eqn{(u[,1]'*x_i*v_i[,1])^2} measures the relative link between
\eqn{R^p} and \eqn{R^{q_i}} associated to xi. It corresponds to a partial squared
singular value notion, since \eqn{\sum_i (u[,1]'*x_i*v_i[,1])^2=s^2},
where s is the usual first singular value of x.
The second solution is obtained from the same criterion, but after
replacing each xi by \code{xi-xi*vi[,1]*vi[,1]^prime}.  And so on for the
successive solutions 1,2,...,r .  The biggest number of solutions may
be r=inf(p,qi), when the xi's are supposed with full rank; then
\code{rmax=min([min(H),p])}.
}
\examples{
x <- matrix(runif(200),10,20)
s <- svdcp(x,c(5,5,10),1)
ss <- svd(x);ss$d[1]^2
sum(s$s2)

}
\references{
Lafosse R. & Hanafi M.(1997) Concordance d'un tableau avec K tableaux: Definition de K+1 uples synthetiques. Revue de Statistique Appliquee vol.45,n.4.
}
\author{
\enc{Lafosse, R.}{R. Lafosse}
}
