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
}
\value{
zero-based vector corresponding to nodes
}
\references{Boost Graph Library  }
\author{ interface by Vince Carey <stvjc@channing.harvard.edu>}

\seealso{strongConnComp  }

\examples{
if (length(agrep("solaris", version[["platform"]]))==0)
{
km <- fromGXL(file(system.file("GXL/kmstEx.gxl",package="graph")))
km@nodes <- c(km@nodes,"F","G","H")
km@edgeL$F <- list(edges=character(0))
km@edgeL$G <- list(edges=8)
km@edgeL$H <- list(edges=7)
km@edgemode <- "undirected"
connectedComp(km)
} else {
 cat("this function is not working in solaris at present,\n")
 cat("please use a windows or linux platform")
 }
}
\keyword{ models }
