\name{match.pairC}
\alias{match.pairC}

\title{Protein molecule matching}
\description{
     Performs pairwise and multipair matching of protein molecules, using the size and shape likelihood for obtaining the matched mean and the Hungarian algorithm to determine the optimal matching.
     }

\usage{
match.pairC(data, SP = 10, n.cores = 8, volume = NULL,
                        jumps = 0, PAM = NULL, gap_open = 0, gap_ext = 0,
                        restrict = TRUE)
}
\arguments{
\item{data}{A list of length at least 2 in the format created by readpdb. Each element is a data frame of dimensions k x 4, where the first column is the amino acid sequence and the rest the 3-dimensional coordinates.}
\item{SP}{If it is a number then automatic selection of starting points take place, else a matrix p x n with p>4.}
\item{n.cores}{The number of cores to be used for parallel computing.}
\item{volume}{If volume = NULL the volume is calculated from the data.}
\item{jumps}{Number of random matching jumps to be performed, if jumps = 0 algorithm runs without random jumps.}
\item{PAM}{PAM matrix to be used for sequence structure matching, if PAM = NULL no PAM matrix is used.}
\item{gap_open}{Gap opening penalty, if gap_open = 0 no gap opening penalty is used.}
\item{gap_ext}{Gap extension penalty, if gap_ext = 0 no gap extension penalty is used.}
\item{restirct}{If TRUE then search for optimal matching happens only on pairs that preserve sequence order.}
}

\value{A list with components
\item{align}{p x m matrix with the matching correspondence of each landmark.}
\item{M}{number of matched landmarks.}
\item{rmsd}{RMSD of the final alignment.}
\item{TMScore}{TMscore of the final alignment.}
\item{SO}{Structure overlap of the final alignment.}
\item{lik}{log - likelihood of the final alignment.}
\item{time}{Total time needed for the alignment in seconds.}
}     
\author{Thodoris Gkolias}
\examples{
data1 = readpdb('2gb1')
data2 = readpdb('1ubq')
data = list(data1,data2)

ali1 = match.pairC(data, SP = 5, restrict = FALSE)
ali2 = match.pairC(data, SP = 7, PAM = PAM[120], gap_open = 2, gap_ext = 10,
          restrict = FALSE)
}
