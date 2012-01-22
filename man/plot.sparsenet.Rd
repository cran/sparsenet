\name{plot.sparsenet}
\alias{plot.sparsenet}
\title{plot coefficients from a "sparsenet" object}
\description{
Produces a series of coefficient profile plots of the coefficient paths for a
fitted \code{"sparsenet"} object.
}
\usage{
\method{plot}{sparsenet}(x, xvar = c("rsq","lambda","norm"), which.gamma=NULL, label = FALSE,...)
}
\arguments{
  \item{x}{fitted \code{"sparsenet"} model}
  \item{xvar}{What is on the X-axis. \code{"rsq"} plots against the
    percent variance explained on the training data,  \code{"lambda"} against the log-lambda
    sequence, and \code{"norm"} plots against the
    L1-norm of the coefficients}
  \item{which.gamma}{sequence numbers of \code{gamma} values to be used
    in the plots; default is all used in the fit}
  \item{label}{If \code{TRUE}, label the curves with variable sequence numbers.}
  \item{\dots}{Other graphical parameters to plot}
}
\details{
A series of coefficient profile plots is produced, one for each
\code{gamma} specified. Users should set up the appropriate layout.
}
\references{
http://www.stanford.edu/~hastie/Papers/Sparsenet/jasa_MFH_final.pdf
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
fit=sparsenet(x,y)
par(mfrow=c(3,3))
plot(fit)
}
\keyword{ lasso }
\keyword{ subset }
