\name{connectedComp}
\alias{connectedComp}
\title{Identify Connected Components in
an Undirected Graph  }
\description{The connected components in
an undirected graph are identified. If the graph is directed then the
weakly connected components are identified.
}
\usage{
connectedComp(g)
}
\arguments{
  \item{g}{ g: graph with \code{edgemode} \dQuote{undirected}}
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
\references{Boost Graph Library}

\author{ Vince Carey <stvjc@channing.harvard.edu>}
\seealso{\code{\link[graph]{connComp}},\code{\link{strongComp}},
  \code{\link[graph]{ugraph}}}

\examples{
km <- fromGXL(file(system.file("GXL/kmstEx.gxl",package="graph")))
km@nodes <- c(km@nodes,"F","G","H")
km@edgeL$F <- list(edges=numeric(0))
km@edgeL$G <- list(edges=8)
km@edgeL$H <- list(edges=7)
ukm <- ugraph(km)
ukm
edges(ukm)
connectedComp(ukm)
}
\keyword{models }
