\name{sparsenet}
\alias{sparsenet}
\title{
Fit a linear model regularized by the nonconvex MC+ sparsity penalty
}
\description{
  Sparsenet uses coordinate descent on the MC+ nonconvex penalty family,
  and fits a surface of solutions over the two-dimensional parameter
  space. This penalty family is indexed by an overall strength paramter lambda
  (like lasso), and a convexity parameter gamma. Gamma = infinity
  corresponds to the lasso, and gamma = 1 best subset.
}
\usage{
sparsenet(x, y, weights, exclude, dfmax = nvars + 1, pmax = min(dfmax *2, nvars),
ngamma = 9, nlambda = 50, max.gamma = 150, min.gamma = 1.000001,
lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04), lambda = NULL,
gamma = NULL, parms = NULL, warm = c("lambda", "gamma", "both"),
thresh = 1e-05, maxit = 1e+06)
}
\arguments{
  \item{x}{Input matrix of nobs x nvars predictors
}
  \item{y}{response vector
}
  \item{weights}{Observation weights; default 1 for each observation
}
  \item{exclude}{Indices of variables to be excluded from the
    model. Default is none.
}
\item{dfmax}{Limit the maximum number of variables in the
    model. Useful for very large \code{nvars}, if a partial path is desired.}
  \item{pmax}{Limit the maximum number of variables ever to be nonzero}
  \item{ngamma}{
Number of gamma values, if \code{gamma} not supplied; default is 9.
}
  \item{nlambda}{ Number of lambda values, if \code{lambda} not
    supplied; default is 50
}
  \item{max.gamma}{Largest gamma value to be used, apart from infinity
    (lasso), if \code{gamma} not supplied; default is 150
}
  \item{min.gamma}{Smallest value of gamma to use, and should be >1;
    default is 1.000001
}
\item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero). The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.0001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.01}.
  A very small value of
  \code{lambda.min.ratio} will lead to a saturated fit in the \code{nobs <
    nvars} case.
}
  \item{lambda}{A user supplied \code{lambda} sequence, in decreasing order. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
    \code{lambda} overrides this. WARNING: use with care. Do not supply
  a single value for \code{lambda} (for predictions after CV use \code{predict()}
  instead).  Supply instead
    a decreasing sequence of \code{lambda} values. \code{sparsenet} relies
  on its warms starts for speed, and its often faster to fit a whole
  path than compute a single fit.}

  \item{gamma}{
Sparsity parameter vector, with 1<gamma<infty. Gamma=1 corresponds to
best-subset regression, gamma=infty to the lasso. Should be given in
decreasing order.
}
  \item{parms}{
    An optional three-dimensional array: 2x ngamma x nlambda.
    Here the user can supply exactly the gamma, lambda pairs that are to
    be traversed by the coordinate descent algorithm. 
}
  \item{warm}{
How to traverse the grid. Default is "lambda", meaning warm starts from
the previous lambda with the same gamma. "gamma" means the opposite,
previous gamma for the same lambda. "both" tries both warm starts, and
uses the one that improves the criterion the most. 
}
\item{thresh}{Convergence threshold for coordinate descent. Each
  coordinate-descent loop continues until the maximum change in the
  objective after any coefficient update is less than \code{thresh}
  times the null Rss. Defaults value is \code{1E-5}.}
    \item{maxit}{Maximum number of passes over the data for all lambda/gamma
  values; default is 10^6.}
}
\details{
  This algorithm operates like \code{glmnet}, with its alpha parameter
    which moves the penalty between lasso and ridge; here gamma moves it
    between lasso and best subset.
  The algorithm traverses the two dimensional gamma/lambda array in a nested loop, with
  decreasing gamma in the outer loop, and decreasing lambda in the inner
  loop. Because of the nature of the MC+ penalty, each coordinate update
  is a convex problem, with a simple two-threshold shrinking scheme:
  beta< lambda set to zero; beta > lambda*gamma leave alone; beta
  inbetween, shrink proportionally. Note that this algorithm ALWAYS
  standardizes the columns of x and y to have mean zero and variance 1
    (using the 1/N averaging) before it computes its fit. The
    coefficients reflect the original scale.
  }
\value{
An object of class \code{"sparsenet"}, with a number of
    components. Mostly one will access the components via generic
    functions
    like \code{coef()}, \code{plot()}, \code{predict()} etc.
    \item{call}{the call that produced this object}
    \item{rsq}{The percentage variance explained on the training data;
      an ngamma x nlambda matrix.}
  \item{jerr}{error flag, for warnings and errors (largely for internal debugging).}
  \item{coefficients}{A coefficient list with ngamma elements; each of
    these is a coefficient list with various components: the matrix beta
    of coefficients, its dimension dim, the vector of intercepts, the lambda sequence, the gamma value,  the sequence
    of df (nonzero coefficients) for each solution.}
  \item{parms}{Irrespective how the parameters were input, the three-way
    array of what was used.}
  \item{gamma}{The gamma values used}
  \item{lambda}{The lambda values used}
  \item{max.lambda}{The entry value for lambda}
}
\references{
http://www.stanford.edu/~hastie/Papers/Sparsenet/jasa_MFH_final.pdf
}
\author{
Rahul Mazumder, Jerome Friedman and Trevor Hastie

Maintainer: Trevor Hastie <hastie@stanford.edu>
}
\seealso{
\code{glmnet} package,  \code{predict}, \code{coef}, \code{print} and \code{plot} methods, and the \code{cv.sparsenet} function.}
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
