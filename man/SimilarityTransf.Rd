\name{SimilarityTransf}
\alias{SimilarityTransf}

\title{Computes  similarity transformations.}

\description{
     Computes the similarity transformations between two matrices.
}

\usage{
     SimilarityTransf(X, Xrot)
}
\arguments{
     \item{X}{A k x m x n array with the original data, k is the number of landmarks, m is the number of dimensions and n is the sample size.}
     \item{Xrot}{A k x m x n array with the rotated data, k is the number of landmarks, m is the number of dimensions and n is the sample size.}
}

\value{A list with components
     \item{location}{The location parameters.}
     \item{R}{The rotation matrices between the original and rotated data.}
}

\author{Thodoris Gkolias}
\examples{
     X <- array(rnorm(20*3*5,10,2), dim = c(20,3,5))
     Xrot <- GpaC(X)$rot
     
     SimilarityTransf(X, Xrot)
}
