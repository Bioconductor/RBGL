\name{mstree.kruskal}
\alias{mstree.kruskal}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{ Kruskal's minimum spanning tree in boost }

\description{ compute the minimum spanning tree (MST) for
a graph and return a representation in matrices }

\usage{
mstree.kruskal(x) 
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{ instance of class graph }
}

\details{
calls to kruskal minimum spanning tree algorithm of Boost graph library
}

\value{
a list \item{edgeList}{a matrix m of dimension 2 by number of edges in the
MST, with m[i,j] the jth node in edge i}
\item{weights}{a vector of edge weights corresponding to the
columns of \code{edgeList}}
\item{nodes}{the vector of nodes of the input graph \code{x}}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ VJ Carey <stvjc@channing.harvard.edu> }


\examples{
con1 <- file(system.file("XML/kmstEx.gxl",package="RBGL"), open="r")
km <- fromGXL(con1)
close(con1)

mstree.kruskal(km)
edgeData(km, "B", "D", "weight") <- 1.1
edgeData(km, "B", "E", "weight") <- .95
mstree.kruskal(km)

con2 <- file(system.file("XML/telenet.gxl",package="RBGL"), open="r")
km2 <- fromGXL(con2)
close(con2)

m <- mstree.kruskal(km2)
print(sum(m[[2]]))
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
