\name{johnson.all.pairs.sp}
\alias{johnson.all.pairs.sp}

\title{compute shortest path distance matrix for all pairs of nodes}

\description{ compute shortest path distance matrix for all pairs of nodes }

\usage{
johnson.all.pairs.sp(g)
}

\arguments{
  \item{g}{ graph object for which edgeMatrix and edgeWeights are defined}
}

\details{
Uses BGL algorithm.
}

\value{
matrix of shortest path lengths, read from row node to col node
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{  }

\seealso{\code{\link{bellman.ford.sp}}, \code{\link{dag.sp}}, \code{\link{dijkstra.sp}}, \code{\link{sp.between}} }

\examples{
z <- fromGXL(file(system.file("dot/joh.gxl", package="RBGL"), open="r"))
johnson.all.pairs.sp(z)
}
\keyword{ models }
