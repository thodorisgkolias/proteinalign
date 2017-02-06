\name{ReadPDB}
\alias{ReadPDB}


\title{Read PDB data from a directory}
\description{
     Reads pdb format data from a local directory
     }

\usage{
ReadPDB(dir, atom = 'CA')
}
\arguments{
\item{dir}{Local directory of the data}
\item{atom}{The atom of the protein chain to be selected.}
}

\value{A matrix with 4 columns. The first is the amino acid for each residue and the other three the 3-dimensional coordinates of the protein molecule.}
    
\author{Thodoris Gkolias}
\examples{
\dontrun{
## Specify the local directory of a .pdb or .atm file
file_dir <- '1a87.pdb'
ReadPDB(dir = file_dir)
}
}