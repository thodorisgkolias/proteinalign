\name{ProteinAlign}
\alias{ProteinAlign}

\title{Protein molecule alignment}
\description{
     Performs pairwise and multipair alignment of protein molecules, using the size and shape likelihood for obtaining the matched mean and the Hungarian algorithm to determine the optimal matching.
     }

\usage{
ProteinAlign(data, SP = 10, n.cores = 8, volume = NULL,
                        jumps = 0, PAM = NULL, gap_open = 0, gap_ext = 0,
                        restrict = FALSE)
}
\arguments{
\item{data}{A list of length at least 2 in the format created by LoadPDB or ReadPDB.}
\item{SP}{Number or matrix, if a number automatic selection of starting points, else a matrix p x n with p>4.}
\item{n.cores}{The number of cores to be used for parallel computing.}
\item{volume}{If volume = NULL the volume is calculated from the data.}
\item{jumps}{Number of random matching jumps to be performed, if jumps = 0 algorithm runs without random jumps.}
\item{PAM}{PAM matrix to be used for sequence structure matching, if PAM = NULL no PAM matrix is used.}
\item{gap_open}{Gap opening penalty, if gap_open = 0 no gap opening penalty is used.}
\item{gap_ext}{Gap extension penalty, if gap_ext = 0 no gap extension penalty is used.}
\item{restrict}{If TRUE then search for optimal matching happens only on pairs that preserve sequence order.}
}

\value{A list with components
\item{align}{p x m matrix with the matching correspondence of each landmark.}
\item{M}{Number of matched landmarks.}
\item{rmsd}{RMSD of the final alignment.}
\item{TMScore}{TMscore of the final alignment.}
\item{SO}{Structure overlap of the final alignment.}
\item{lik}{log - likelihood of the final alignment.}
\item{time}{Total time needed for the alignment in seconds.}
}     
\author{Thodoris Gkolias}
\examples{
## Loading data from the RCSB protein data bank
data1 <- LoadPDB('2gb1')
data2 <- LoadPDB('1ubq')
data <- list(data1,data2)

##  Alignment using 5 automatic selected landmarks 
ali1 <- ProteinAlign(data, SP = 5, n.cores = 2)

## Alignment using user specified starting points and PAM matrix
SP_user <- cbind(c(2,3,4,5,16,18,22,50), c(2,3,4,5,13,15,19,65))
PAM120 <- PAM[[120]]
ali2 <- ProteinAlign(data, SP = SP_user, n.cores = 2, PAM = PAM120)
}
