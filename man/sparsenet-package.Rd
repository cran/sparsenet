\name{sparsenet-package}
\alias{sparsenet-package}
\docType{package}
\title{
Fit a linear model regularized by the nonconvex MC+ sparsity penalty
}
\description{
  Sparsenet uses coordinate descent on the MC+ nonconvex penalty family,
  and fits a surface of solutions over the two-dimensional parameter space.
}
\details{
At its simplest, provide \code{x,y} data and it returns the solution
paths. There are tools for prediction, cross-validation, plotting and printing.

}
\author{     
Rahul Mazumder, Jerome Friedman and Trevor Hastie

Maintainer: Trevor Hastie <hastie@stanford.edu>
}
\references{
  Mazumder, Rahul, Friedman, Jerome and Hastie, Trevor (2011)
  \emph{SparseNet: Coordinate Descent with Nonconvex Penalties. JASA, Vol 106(495), 1125-38},
 \url{https://hastie.su.domains/public/Papers/Sparsenet/Mazumder-SparseNetCoordinateDescent-2011.pdf}
}
\keyword{ regression }
\keyword{ package }
\keyword{ sparse }
\examples{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
fit=sparsenet(x,y)
plot(fit)
cvfit=cv.sparsenet(x,y)
plot(cvfit)
}
