\name{RMSD}
\alias{RMSD}

\title{Root Mean Square Deviation}
\description{
     Computes the Root Mean Square Deviation between two matrices.
}

\usage{
     RMSD(X, Y)
     
}
\arguments{
     \item{X}{A k x m x n matrix.}
     \item{Y}{A k x m x n matrix.}
}

\value{The RMSD value between X and Y.
}

\author{Thodoris Gkolias}
\examples{
     ## Create data
     X <- matrix(rnorm(30), 10, 3)
     Y <- matrix(rnorm(30, 3, 2), 10, 3)
     RMSD(X, Y)
}
