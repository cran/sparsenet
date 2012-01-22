\name{gendata}
\alias{gendata}
\title{
Generate data for testing sparse model selection
}
\description{
This function generates x/y data for testing sparsenet and glmnet
}
\usage{
gendata(N, p, nonzero, rho, snr = 3, alternate = TRUE)
}
\arguments{
  \item{N}{
Sample size (eg 500)
}
  \item{p}{
Number of features or variables (eg 1000)
}
  \item{nonzero}{
Number if nonzero coefficients (eg 30)
}
  \item{rho}{
pairwise correlation between features
}
  \item{snr}{
Signal to noise ratio - SD signal/ SD noise - try 3
}
  \item{alternate}{
Alternate sign of coefficients
}
}
\details{
Generates Gaussian x and y data. The nonzero coefficients decrease
linearly in absolute value from nonzero down to 0. If
\code{alternate=TRUE} their signs alternate, else not
}
\value{
A list with components x and y as well some other details about the dataset
}
\author{
Trevor Hastie and Jerome Friedman
}
\examples{
train.data=gendata(100,1000,nonzero=30,rho=0.3,snr=3)
fit=sparsenet(train.data$x,train.data$y)
par(mfrow=c(3,3))
plot(fit)
par(mfrow=c(1,1))
fitcv=cv.sparsenet(train.data$x,train.data$y,trace.it=TRUE)
plot(fitcv)
}
\keyword{ simulate }
\keyword{ regression}
\keyword{ subset }
