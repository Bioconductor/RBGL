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
}
\value{
list with element 1 the connectivity, and element 2 a list
of pairs of nodes defining the edges of the minimum
disconnecting set
}
\references{Boost Graph Library  }
\author{ interface by Vince Carey <stvjc@channing.harvard.edu>}

\seealso{connectedComp  }

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
edgeConnectivity(coex)
coex2 <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
edgeConnectivity(coex2)
}
\keyword{ models }
