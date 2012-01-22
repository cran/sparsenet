\name{predict.cv.sparsenet}
\alias{coef.cv.sparsenet}
\alias{predict.cv.sparsenet}
\title{make predictions from a "cv.sparsenet" object.}
\description{
This function makes predictions from a cross-validated sparsenet model,
using the stored \code{"sparsenet.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.sparsenet}(object, newx, which=c("parms.min","parms.1se"),...)
\method{coef}{cv.sparsenet}(object, which=c("parms.min","parms.1se"),...)
}
\arguments{
  \item{object}{Fitted \code{"cv.sparsenet"}  object.}
  \item{newx}{Matrix of new values for \code{x} at which predictions are
  to be made. Must be a matrix. See documentation for \code{predict.sparsenet}.}
  \item{which}{Either the paramaters of the minimum of the CV curves
    (default \code{"parms.min"} or the parameters corresponding to the
    one standard-error rule \code{parms.1se})}
  \item{\dots}{Not used. Other arguments to predict. }
}
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on
to the \code{predict} method for \code{sparsenet} objects.}
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
fitcv=cv.sparsenet(x,y)
predict(fitcv,x)
}
\keyword{ lasso }
\keyword{ subset }
