\name{connectedComp}
\alias{connectedComp}
\title{identify connected components in
an undirected graph  }
\description{identify connected components in
an undirected graph
}
\usage{
connectedComp(g)
}
\arguments{
  \item{g}{ g: graph with edgemode undirected}
}
\details{
Uses a DFS based approach.
Mathematical formalism can be
found in graph package
\code{\link{connComp}}, see also \code{\link{strongComp}}.  
}

\value{
List of node names defining the connected components.
}
\references{Boost Graph Library  }

\author{ interface by Vince Carey <stvjc@channing.harvard.edu>}

\examples{
km <- fromGXL(file(system.file("GXL/kmstEx.gxl",package="graph")))
km@nodes <- c(km@nodes,"F","G","H")
km@edgeL$F <- list(edges=character(0))
km@edgeL$G <- list(edges=8)
km@edgeL$H <- list(edges=7)
ukm <- ugraph(km)
ukm
edges(ukm)
connectedComp(ukm)
}
\keyword{models }
