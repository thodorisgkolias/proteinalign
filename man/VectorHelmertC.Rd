\name{VectorHelmertC}
\alias{VectorHelmertC}

\title{Helmertized landmarks}
\description{
     Produces a matrix with the Helmertized landmarks of X 
}

\usage{
     VectorHelmertC(X)
     
}
\arguments{
     \item{X}{A k x m x n matrix with the data.}
}

\value{A (k-1) x m x n with the Helmertized landmarks of X.
}

\author{Alfred Kume}
\examples{
     ## Create data
     x <- array(rnorm(3*10*5), dim=c(10, 3, 5))
     
     VectorHelmertC(x)
}
