\name{cv.sparsenet}
\alias{cv.sparsenet}
\title{Cross-validation for sparsenet}
\description{Does k-fold cross-validation for sparsenet, produces a plot,
  and returns values for \code{gamma, lambda}}
\usage{
cv.sparsenet(x, y, weights, type.measure = c("mse", "mae"), ...,nfolds = 10,foldid, trace.it=FALSE)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{sparsenet}.}
  \item{y}{response \code{y} as in \code{sparsenet}.}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{type.measure}{loss to use for cross-validation. Currently two
    options:
    squared-error (\code{type.measure="mse"}) or
mean-absolute error ( \code{type.measure="mae"} )
 } 
  \item{\dots}{Other arguments that can be passed to \code{sparsenet}.}
  \item{nfolds}{number of folds - default is 10. Although \code{nfolds}
    can be as large as the sample size (leave-one-out CV), it is not
    recommended for large datasets. Smallest value allowable is \code{nfolds=3}}
  \item{foldid}{an optional vector of values between 1 and \code{nfold}
    identifying whhat fold each observation is in. If supplied,
    \code{nfold} can be missing.}
  \item{trace.it}{If \code{TRUE}, then we get a printout that shows the
    progress}
  }
\details{The function runs \code{sparsenet} \code{nfolds}+1 times; the
  first to get the \code{lambda} sequence, and then the remainder to
  compute the fit with each of the folds omitted. The error is
  accumulated, and the average error and standard deviation over the
  folds is computed. 
}
\value{an object of class \code{"cv.sparsenet"} is returned, which is a
  list with the ingredients of the cross-validation fit.
  \item{lambda}{the values of \code{lambda} used in the fits. This is an
  \code{nlambda x ngamma} matrix}
  \item{cvm}{The mean cross-validated error - a matrix shaped like lambda}
  \item{cvsd}{estimate of standard error of \code{cvm}.}
  \item{cvup}{upper curve = \code{cvm+cvsd}.}
  \item{cvlo}{lower curve = \code{cvm-cvsd}.}
  \item{nzero}{number of non-zero coefficients at each \code{lambda,
      gamma} pair.}
  \item{name}{a text string indicating type of measure (for plotting
    purposes).}
  \item{sparsenet.fit}{a fitted sparsenet object for the full data.}
  \item{call}{The call that produced this object}
  \item{parms.min}{values of \code{gamma, lambda} that gives minimum
    \code{cvm}.}
  \item{which.min}{indices for the above}
  \item{lambda.1se}{\code{gamma, lambda} of smallest model (df) such that error is
    within 1 standard error of the minimum.}
  \item{which.1se}{indices of the above}
}
\references{
http://www.stanford.edu/~hastie/Papers/Sparsenet/jasa_MFH_final.pdf
}
\author{
Rahul Mazumder, Jerome Friedman and Trevor Hastie

Maintainer: Trevor Hastie <hastie@stanford.edu>
}
\seealso{
\code{glmnet} package,  \code{predict}, \code{coef}, \code{print} and \code{plot} methods, and the \code{sparsenet} function.}
\examples{
train.data=gendata(100,1000,nonzero=30,rho=0.3,snr=3)
fit=sparsenet(train.data$x,train.data$y)
par(mfrow=c(3,3))
plot(fit)
par(mfrow=c(1,1))
fitcv=cv.sparsenet(train.data$x,train.data$y,trace.it=TRUE)
plot(fitcv)
}
\keyword{ lasso }
\keyword{ subset }
