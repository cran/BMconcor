% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concorsreg.R
\name{concorsreg}
\alias{concorsreg}
\title{Redundancy of sets yj by one set x}
\usage{
concorsreg(x, px, y, py, r)
}
\arguments{
\item{x}{are the \code{n} times \code{p} and \code{n} times \code{q} matrices of \code{p} and \code{q} centered column}

\item{px}{The row vector which contains the numbers pi, i = 1,...,kx, of the kx subsets xi of x : \eqn{\sum_i p_i}=sum(px)=p. px is the partition vector of x}

\item{y}{See \code{x}}

\item{py}{The partition vector of y. A row vector containing the numbers \code{qi} for \code{i = 1,...,ky} of the \code{ky} subsets \code{yi} of \code{y : sum(qi)=sum(py)=q}.}

\item{r}{The number of wanted successive solutions}
}
\value{
A \code{list} with following components:
\item{cx}{a \code{n} times \code{r}matrix of the r explanatory components}
\item{v}{is a \eqn{q \times r}{q x r} matrix of ky row blocks \eqn{v_i} (\eqn{q_i \times r}{q_i x r}) of axes in Rqi relative to yi; \eqn{v_i'*v_i = \mbox{Id}}}
\item{varexp}{is a \eqn{ky \times r}{ky x r} matrix; each column k contains ky explained variances \eqn{\rho(cx[,k],y_i*v_i[,k])^2 \mbox{var}(y_i*v_i[,k])}}
}
\description{
Regression of several subsets of variables Yj by another set X. SUCCESSIVE SOLUTIONS
}
\examples{

x <- matrix(runif(50),10,5);y <- matrix(runif(90),10,9)
x <- scale(x);y <- scale(y)
crs <- concorsreg(x,c(2,3),y,c(3,2,4),2)
crs$varexp[1,1,]

}
\author{
\enc{Lafosse, R.}{R. Lafosse}
}
