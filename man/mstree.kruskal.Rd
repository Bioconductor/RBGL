\name{mstreeKruskal}
\alias{mstreeKruskal}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Kruskal's minimum spanning tree in boost }
\description{
compute the minimum spanning tree (MST) for
a graph and return a representation in matrices
}
\usage{
mstreeKruskal(x) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ instance of class graph }
}
\details{
calls to kruskal minimum spanning tree algorithm of Boost graph library
}
\value{
a list
\item{edgeList}{a matrix m of dimension 2 by number of edges in the
MST, with m[i,j] the jth node in edge i}
\item{weights}{a vector of edge weights corresponding to the
columns of \code{edgeList}}
\item{nodes}{the vector of nodes of the input graph \code{x}}
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\examples{
km <- fromGXL(file(system.file("GXL/kmstEx.gxl",package="graph")))
mstreeKruskal(km)
km@edgeL$B$weights <- c(1.1,.95)
mstreeKruskal(km)
km2 <- fromGXL(file(system.file("XML/telenet.gxl",package="RBGL")))
m <- mstreeKruskal(km2)
print(sum(m[[2]]))
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
