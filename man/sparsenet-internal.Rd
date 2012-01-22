\name{sparsenet-internal}
\title{Internal sparsenet functions}
\alias{argmin}
\alias{getcoef_list}
\alias{lambda0}
\alias{sparsepredict}
\alias{zeromat}
\alias{summary.sparsenet}
\alias{print.cv.sparsenet}
\alias{print.sparsenet}
\description{Internal sparsenet functions}
\author{Trevor Hastie}
\details{These are not intended for use by users. \code{argmin}
  identifies the lambda, gamma pair that minimizes cv
  error. \code{lambda0} should go away; it currently evaluates the
  entry value for lambda.
}
\keyword{internal}
