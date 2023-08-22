\name{plot.cv.sparsenet}
\alias{plot.cv.sparsenet}
\title{plot the cross-validation curves produced by cv.sparsenet}
\description{
Plots the cross-validation curves for each value of \code{gamma} in one figure, as a function of the \code{lambda} values used.}
\usage{
\method{plot}{cv.sparsenet}(x,  ...)
}
\arguments{
  \item{x}{fitted \code{"cv.sparsenet"} object}
  \item{\dots}{Other graphical parameters to plot}
}
\details{
A plot is produced, and nothing is returned.}
\references{
 Mazumder, Rahul, Friedman, Jerome and Hastie, Trevor (2011)
  \emph{SparseNet: Coordinate Descent with Nonconvex Penalties. JASA, Vol 106(495), 1125-38},
 \url{https://hastie.su.domains/public/Papers/Sparsenet/Mazumder-SparseNetCoordinateDescent-2011.pdf}
}


\author{
Rahul Mazumder, Jerome Friedman and Trevor Hastie

Maintainer: Trevor Hastie <hastie@stanford.edu>
}
\seealso{
\code{glmnet} package, \code{sparsenet},  \code{cv.sparsenet}  and
\code{print} and \code{plot} methods for both.}
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fitcv=cv.sparsenet(x,y)
plot(fitcv)
}
\keyword{ lasso }
\keyword{ subset }
