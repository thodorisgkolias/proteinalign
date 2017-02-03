\name{Norm}
\alias{Norm}


\title{Euclidean norm}
\description{
     Computes the Euclidean norm of a matrix 
}

\usage{
     Norm(X)
}
\arguments{
     \item{X}{A matrix k x m}
}

\value{The Euclidean norm of X}

\author{Alfred Kume}
\examples{
x = rnorm(15,5,3)
Norm(x)
}
