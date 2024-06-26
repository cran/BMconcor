% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svdbip.R
\name{svdbip}
\alias{svdbip}
\title{SVD for one bipartitioned matrix x}
\usage{
svdbip(x, K, H, r)
}
\arguments{
\item{x}{a \code{p} times \code{q} matrix}

\item{K}{is a row vector which contains the numbers pk, k=1,...,kx, of the partition of x with kx row blocks : \code{sum(pk)=p}}

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
SVD for bipartitioned matrix x. r successive Solutions
}
\details{
The first solution calculates kx+ky normed vectors: kx vectors \code{uk[:,1]}
of \eqn{R^{p_k}} associated to ky vectors \code{vh[:,1]}'s of \eqn{R^{q_h}}, by
maximizing \eqn{\sum_k \sum_h (u_k[:,1]^prime*x_{kh}*v_h[:,1])^2}, with kx+ky
norm constraints. A value \eqn{(u_k[,1]^prime*x_{kh}*v_h[,1])^2} measures the
relative link between \eqn{R^{p_k}} and \eqn{R^{q_h}} associated to the block xkh.
The second solution is obtained from the same criterion, but after
replacing each xhk by xkh-xkh\emph{vh}vh'-uk\emph{uk'xkh+uk}uk'xkh\emph{vh}vh'.  And
so on for the successive solutions 1,2,...,r .  The biggest number of
solutions may be r=inf(pk,qh), when the xkh's are supposed with full
rank; then \code{rmax=min([min(K),min(H)])}.
When K=p (or H=q, with t(x)), svdcp function is better.  When H=q and
K=p, it is the usual svd (with squared singular values).
Convergence of algorithm may be not global. So the below proposed
initialisation of the algorithm may be not very suitable for some data
sets.  Several different random initialisations with normed vectors
might be considered and the best result then choosen.
}
\examples{

x <- matrix(runif(200),10,20)
s <- svdbip(x,c(3,4,3),c(5,15),3)

}
\references{
Kissita G., Cazes P., Hanafi M. & Lafosse (2004) Deux methodes d'analyse factorielle du lien entre deux tableaux de variables partitiones. Revue de Statistique Appliquee.
}
\author{
\enc{Lafosse, R.}{R. Lafosse}
}
