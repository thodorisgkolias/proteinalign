\name{SpChoose}
\alias{SpChoose}

\title{Starting point selection.}

\description{
     Selects a starting alignment between two protein molecules, using sequence and structural alignment.
}

\usage{
     SpChoose(data1, data2, k = 4, cut = 0.2)

}
\arguments{
     \item{data1}{A k x 4 matrix with amino acid sequence information in the fist column and the 3-dimensional coordinates in the rest three.}
     \item{data2}{A l x 4 matrix with amino acid sequence information in the fist column and the 3-dimensional coordinates in the rest three.}
     \item{k}{Number of starting points to be selected.}
     \item{cut}{Cutoff TMscore under which the starting points are selected.}
}

\value{A k x 2 matrix with the selected starting points for the alignment.
}

\author{Thodoris Gkolias}
\examples{
     ## Load data
     data1 <- LoadPDB('1acx')
     data2 <- LoadPDB('1cob')
     SpChoose(data1, data2)
}
