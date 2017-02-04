\name{mclapply.win}
\alias{mclapply.win}

\title{mclapply function for Windows}
\description{
     A wrapper for the mclapply function to be used in Windows systems, using the parLapply function.
}

\usage{
    mclapply.win(X,FUN,mc.cores,...)
}
\arguments{
     \item{X}{a vector (atomic or list) or an expressions vector.}
     \item{FUN}{the function to be applied to (parLapply) each element of X.}
     \item{mc.cores}{Number of cores to use.}
     \item{...}{Optional arguments to FUN.}
}

\value{A list of the same length as X and named by X.}

\author{Thodoris Gkolias}
\examples{
    mclapply.win(1:4, function(xx){ Sys.sleep(3) }, mc.cores = 2 )
     
}
