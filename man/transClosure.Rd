\name{transitive.closure}
\alias{transitive.closure}
\title{Compute transitive closure of a graph}
\description{Compute  transitive closure of a graph}
}
\usage{
transitive.closure(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class }
}
\details{
The transitive closure of a graph G = (V,E) is a graph G* = (V,E*) such that E* 
contains an edge (u,v) if and only if G contains a path (of at least one edge) 
from u to v. 
}
\value{
\item{vertices}{the vertices of the transitive closure }
\item{edges}{the edges of the transitive closure }
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
transitive.closure(coex)
}
\keyword{ models }
