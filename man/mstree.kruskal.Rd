\name{mstree.kruskal}
\alias{mstree.kruskal}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Kruskal's minimum spanning tree in boost }
\description{
returns list of two matrices: 2XNE of edges of MST and
1XNE of weights of these edges, and of the vector of nodes
}
\usage{
mstree.kruskal(x) # now x assumed to be Bioconductor graph graphNEL
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ instance of class graphNEL from Bioconductor graph class}
}
\details{
calls to kruskal minimum spanning tree
}
\value{
list of two matrices and vector of nodes, see above
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\examples{
km <- fromGXL(file(system.file("GXL/kmstEx.gxl",package="graph")))
mstree.kruskal(km)
km@edgeL$B$weights <- c(1.1,.95)
mstree.kruskal(km)
km2 <- fromGXL(file(system.file("XML/telenet.gxl",package="RBGL")))
m <- mstree.kruskal(km2)
print(sum(m[[2]]))
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
