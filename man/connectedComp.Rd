\name{connectedComp}
\alias{connectedComp}
\title{Identify Connected Components in an Undirected Graph  }
\description{The connected components in an undirected graph are identified. 
If the graph is directed then the weakly connected components are identified. }

\usage{
connectedComp(g)
}

\arguments{
  \item{g}{ graph with \code{edgemode} \dQuote{undirected}}
}

\details{
Uses a depth first search approach to identifying all the connected
components of an undirected graph. If the input, \code{g}, is a directed
graph it is first transformed to an undirected graph (using
\code{\link[graph]{ugraph}}).
}

\value{
  A list of length equal to the number of connected components in
  \code{g}. Each element of the list contains a vector of the node
  labels for the nodes that are connected.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ Vince Carey <stvjc@channing.harvard.edu>}

\seealso{\code{\link[graph:graph-class]{connComp}},\code{\link{strongComp}}, \code{\link[graph]{ugraph}}, \code{\link{same.component}} }

\examples{
con <- file(system.file("GXL/kmstEx.gxl",package="graph"), open="r")
km <- fromGXL(con)
close(con)
km <- graph::addNode(c("F","G","H"), km)
km <- addEdge("G", "H", km, 1)
km <- addEdge("H", "G", km, 1)
ukm <- ugraph(km)
ukm
edges(ukm)
connectedComp(ukm)
}
\keyword{models }
