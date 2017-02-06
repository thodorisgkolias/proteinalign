\name{LogLikSizeShape}
\alias{LogLikSizeShape}

\title{Size and shape loglikelihood}
\description{
     Computes the Normal matching size and shape log-likehood for a matrix X inside a volume A.
     }

\usage{
LogLikSizeShape(X, mu, sigma, p, M, A)
}
\arguments{
\item{X}{3 x k matrix with the data}
\item{mu}{3 x k matrix with the mean}
\item{sigma}{standard deviation of the data}
\item{p}{Number of the matched points of X}
\item{M}{Number of landmarks of X}
\item{A}{Volume where the data are contained}
}

\value{The value of the size and shape log-likelihood}
  
\author{Thodoris Gkolias}
\examples{
x <- matrix(rnorm(30), 3, 10)
meanmat <- matrix(rnorm(30,2,1), 3, 10)

LogLikSizeShape(X = x, mu = meanmat, sigma = 1, p = 5, M = 10, A = 1000)
}