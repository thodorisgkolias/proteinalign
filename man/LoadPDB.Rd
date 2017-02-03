\name{LoadPDB}
\alias{LoadPDB}

\title{Load PDB data}
\description{
     Loads data from the RCSB protein data bank
     }

\usage{
LoadPDB(pdb, atom = 'CA', chain = 'A')
}
\arguments{
\item{pdb}{The name of the protein molecule to be loaded}
\item{atom}{The atom of the protein chain to be selected.}
\item{chain}{The chain of the protein to be selected.}
}

\value{A matrix with 4 columns. The first is the amino acid for each residue and the other three the 3-dimensional coordinates of the protein molecule.}
    
\author{Thodoris Gkolias}
\examples{
## Loading data from the RCSB protein data bank
LoadPDB('1acx')
LoadPDB('1cob', chain = 'B')
}
