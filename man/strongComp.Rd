\name{strongComp}
\alias{strongComp}
\title{identify strongly connected components in
a directed graph  }
\description{identify strongly connected components in
an directed graph
}
\usage{
strongComp(g)
}
\arguments{
  \item{g}{ g: graph with edgemode directed}
}
\details{
Uses Tarjan's algorithm.
}
\value{
zero-based vector corresponding to nodes
}
\references{Boost Graph Library  }
\author{ interface by Vince Carey <stvjc@channing.harvard.edu>}

\seealso{connectedComp  }

\examples{
km <- fromGXL(file(system.file("XML/kmstEx.gxl",package="RBGL")))
km@nodes <- c(km@nodes,"F","G","H")
km@edgeL$F <- list(edges=character(0))
km@edgeL$G <- list(edges=8)
km@edgeL$H <- list(edges=7)
strongComp(km)
connectedComp(ugraph(km))
}
\keyword{ models }
