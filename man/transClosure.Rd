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
This function calculates the transitive closure of a given graph.  
See documentation on this function in Boost Graph Library for more details.
}
\value{
  A list of
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
