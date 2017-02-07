\name{RiemdistC}
\alias{RiemdistC}

\title{Riemannian distance}
\description{
     Calculates the riemannian distance for two configurations.
}

\usage{
     RiemdistC(x, y, preshape = 0, reflect = 0)

}
\arguments{
     \item{x}{A k x m matrix.}
     \item{y}{A k x m matrix.}
     \item{preshape}{If 1 then the preshapes of x and y are used.}
     \item{reflect}{If 1 then reflection information is used.}
}

\value{The riemannian distance between x and y
}

\author{Alfred Kume}
\examples{
     ## Create data
     x <- matrix(rnorm(30), 10, 3)
     y <- matrix(rnorm(30,5,1), 10 ,3)
     RiemdistC(x = x, y = y)
}
