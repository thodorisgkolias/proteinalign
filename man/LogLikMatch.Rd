\name{LogLikMatch}
\alias{LogLikMatch}

\title{Matching size and shape log likelihood}
\description{
     Computes the log likelihood of two or more protein molecules being observations from the same mean.
}

\usage{
     LogLikMatch(pair = "NULL", data, seq, matched, vol, PAM, gap_open, gap_ext)
}
\arguments{
     \item{pair}{If not NULL a vector of the same length as the number of protein molecules. Then, the fucntion computes the log likehood of the protein molecules with this pair considered as matched. }
     \item{data}{A list of length at least 2 with the 3-dimensional coordinates of the protein molecules.}
     \item{seq}{A list with the same length as data with the amino acid sequence for each protein molecule.}
     \item{matched}{A p x n matrix, where n is the length of data, with the matched landmarks between the protein molecules.}
     \item{vol}{The volume to be used for the log likelihood.}
     \item{PAM}{PAM matrix to be used for sequence structure matching, if NULL no PAM matrix is used.}
     \item{gap_open}{Gap opening penalty, if 0 no gap opening penalty is used.}
     \item{gap_ext}{Gap extension penalty, if 0 no gap extension penalty is used.}
}

\value{The log likelihood of the protein molecules to come from the same mean.}

\author{Thodoris Gkolias}
\examples{
     ## Load data
     data1 <- LoadPDB('2gb1')
     data2 <- LoadPDB('4fxc')
     data <- list(data1[,-1],data2[,-1])
     seqs <- list(data1[,1],data2[,1])
     
     ## Matched pairs
     matched_pairs <- cbind(1:30,1:30)
     
     ## Pair to be considered as a new match
     new <- c(35,96)
     
     ## Log likelihood of the two matched molecules to come from the same mean.
     LogLikMatch(pair = NULL, data = data, seq = seqs, matched = matched_pairs,
     vol = 30000, PAM = NULL, gap_open = 0, gap_ext = 0)
     
     ## Log likelihood of the two matched molecules with the new pair of landmarks
     ## to come from the same mean.
     LogLikMatch(pair = new, data = data, seq = seqs, matched = matched_pairs,
     vol = 30000, PAM = NULL, gap_open = 0, gap_ext = 0)
}
