\name{EMSizeShape}
\alias{EMSizeShape}

\title{EM algorithm for the size and shape likelihood.}
\description{
     An EM algorithm for estimating the missing rotations of partially observed size and shape data. An isotropic covariance matrix is considered and the helmertized landmarks are used to standardize for location. 
}

\usage{
     EMSizeShape(x)
}
\arguments{
     \item{x}{An array k x m x n where k is the number of landmarks, m the dimensions (2 or 3) and n the number of observations.}
}

\value{A list with components
     \item{mean}{A m x (k-1) with the size and shape mean.}
     \item{sigma}{The standard deviation.}
     \item{lik}{The size and shape log likelihood.}
     \item{rot.data}{A m x (k-1) x n array with the rotated data X}
}

\author{Thodoris Gkolias}
\examples{
     ## Loading data from the RCSB protein data bank
     LoadPDB('1acx')
     LoadPDB('1cob', chain = 'B')
}
