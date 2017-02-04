\name{PAMLogLik}
\alias{PAMLogLik}

\title{Amino acid sequence log likelihood.}
\description{
    Computes the amino acid sequence log likehood between the matched residues for a given PAM matrix.
}

\usage{
     PAMLogLik(seq1, seq2, p, PAM)
}
\arguments{
     \item{seq1}{A k-length vector with the amino acids for the first sequence. The matched amino acids should be first.}
     \item{seq2}{A l-length vector with the amico acids for the second sequence. The matched amino acids should be first.}
     \item{p}{Number of matched amino acids between the two sequences.}
     \item{PAM}{The PAM matrix to be used}
}

\value{The log likelihood for the two matched amino acid sequences.}

\author{Thodoris Gkolias}
\examples{
     ## Amino acid sequences
     aa1 = c('A','A','V','G','M','K','N','N')
     aa2 = c('A','P','G','G','N','P','D','E','E','Q')
     
     PAMLogLik(seq1 = aa1, seq2 = aa2, p = 5, PAM = PAM[[70]])
     
}
