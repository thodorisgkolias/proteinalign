\name{SvdNew}
\alias{SvdNew}


\title{Singular Value Decomposition of a matrix}
\description{
     A wrapper of the svd function, computing the singular value decomposition of a matrix and adjusting for the reflection information.
}

\usage{
     SvdNew(X)
}
\arguments{
     \item{X}{A matrix k x k.}
}

\value{A list which contains
\item{d}{a vector containing the singular values of X}
\item{u}{a matrix whose columns contain the left singular vectors of X}
\item{v}{a matrix whose columns contain the right singular vectors of X}
}

\author{Alfred Kume}
\examples{
     x <- matrix(rnorm(15,5,3))
     Norm(t(x) \%*\% x)
}
