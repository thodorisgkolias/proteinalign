\name{GpaC}
\alias{GpaC}

\title{General Procrustes Analysis in C.}
\description{
     Performs the General Procrustes Analysis written in C.
}

\usage{
     GpaC(x, rescale = 1, reflect = 0, tol1 = 1e-05, tol2 = tol1)
}
\arguments{
     \item{x}{An array k x m x n where k is the number of landmarks, m the dimensions (2 or 3) and n the number of observations.}
     \item{rescale}{If 1 then scaling information is used.}
     \item{reflect}{If 1 then reflection information is used.}
     \item{tol1}{Tolerance for the rotation updates.}
     \item{tol2}{Tolerance for the rescaling updates.}
}

\value{A list with components
     \item{rot}{A k x m x n array with the rotated data of X}
     \item{mshape}{A k x m matrix with the Procrustes mean.}
}

\author{Alfred Kume}
\examples{
     ## Create data
     data <- array(rnorm(20*3*5,mean = 2, sd = 1), dim = c(20,3,5))
     GpaC(data)
}
