\name{edgeConnectivity}
\alias{edgeConnectivity}
\title{computed edge connectivity and min disconnecting set
for an undirected graph}
\description{computed edge connectivity and min disconnecting set
for an undirected graph}
}
\usage{
edgeConnectivity(g)
}
\arguments{
  \item{g}{ g: graph with edgemode undirected}
}
\details{
Consider a graph G consisting of a single connected component.
The edge connectivity of G is the minimum number
of edges in G that can be cut to produce a graph with two
(disconnected) components.  The set of edges in this cut
is called the minimum disconnecting set.
}
\value{
list with element 1 the connectivity, and element 2 a list
of pairs of nodes defining the edges of the minimum
disconnecting set
}
\references{Boost Graph Library by Siek et al.  }
\author{ interface by Vince Carey <stvjc@channing.harvard.edu>}

\seealso{connectedComp  }

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
edgeConnectivity(coex)
coex2 <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
edgeConnectivity(coex2)
}
\keyword{ models }
