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
\tabular{ll}{
Package: \tab sparsenet\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-22\cr
License: \tab GPL-2\cr
}
At its simplest, provide \code{x,y} data and it returns the solution
paths. There are tools for prediction, cross-validation, plotting and printing.

}
\author{     
Rahul Mazumder, Jerome Friedman and Trevor Hastie

Maintainer: Trevor Hastie <hastie@stanford.edu>
}
\references{
http://www.stanford.edu/~hastie/Papers/Sparsenet/jasa_MFH_final.pdf
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
