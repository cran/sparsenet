\name{predict.sparsenet}
\alias{coef.sparsenet}
\alias{predict.sparsenet}
\title{make predictions from a "sparsenet" object.}
\description{
Similar to other predict methods, this functions predicts fitted values,
coefficients and more from a fitted \code{"sparsenet"} object.
}
\usage{
\method{predict}{sparsenet}(object, newx, s = NULL,  which.gamma = NULL,
type=c("response","coefficients","nonzero"), exact = FALSE, ...)
\method{coef}{sparsenet}(object,s=NULL, which.gamma = NULL,exact=FALSE, ...)
}
\arguments{
  \item{object}{Fitted \code{"sparsenet"} model object.}
  \item{newx}{Matrix of new values for \code{x} at which predictions are
  to be made. Must be a matrix. This argument is not used for \code{type=c("coefficients","nonzero")}}
  \item{s}{Value(s) of the penalty parameter \code{lambda} at which
    predictions are required. Default is the entire sequence used to
    create the model. }
  \item{which.gamma}{Index or indices of \code{gamma} values at which
    predictions are to be made. Default is all those used in the fit}
  \item{type}{\code{"response"} returns fitted predictions at
  \code{newx}.  Type
    \code{"coefficients"} computes the coefficients at the requested
    values for \code{s}.   Type  \code{"nonzero"} returns lists of the indices of the nonzero
    coefficients for each value of \code{s}.}
  \item{exact}{By default (\code{exact=FALSE}) the predict function uses linear interpolation
  to make predictions for values of \code{s} that do not coincide with
  those used in the fitting algorithm. Currently \code{exact=TRUE} is
  not implemented, but prints an error message telling the user how to
  achieve the exact predictions. This is done my rerunning the algorithm
with the desired values interspersed (in order) with the values used in
the original fit}
  \item{\dots}{Not used. Other arguments to predict. } 
}
\details{The shape of the objects returned depends on which
  \code{which.gamma} has more than one element.
  If more than one element, a list of predictions is returned, one for
  each gamma.
  }
\value{The object returned depends on type.}
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
predict(fit, which.gamma=5,type="nonzero")
predict(fit,x)
}
\keyword{ lasso }
\keyword{ subset }
