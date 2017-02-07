\name{SimilarityMetrics}
\alias{SimilarityMetrics}

\title{Similarity measures}
\description{
     Computes similarity measures to assess the alignment of two protein molecules.
}

\usage{
     SimilarityMetrics(data1, data2, sol) 
}
\arguments{
\item{data1}{A k x 3 matrix with the 3-dimensional coordinates.}
\item{data2}{A l x 3 matrix with the 3-dimensional coordinates.}
\item{sol}{A p x 2 matrix with alignment between data1 and data2.}
}

\value{A list with components
     \item{rmsd}{The RMSD between the matched parts of data1 and data2 after they have been optimal rotated.}
     \item{TMscore}{Template Modelling score between the mathced parts of data1 and data2 after they have been optimal rotated.}
     \item{SO}{Stucture Overlap between the matched parts of data1 and data2.}
}

\author{Thodoris Gkolias}
\examples{
     ## Load data
     x <- LoadPDB('2gb1')
     x <- as.matrix(x[,-1])
     y <- LoadPDB('1ubq')
     y <- as.matrix(y[,-1])
     ## Alignment between x and y
     ali <- cbind(1:25, 1:25)
     SimilarityMetrics(data1 = x, data2 = y, sol = ali)
}
