\name{floyd.warshall.all.pairs.sp}
\alias{floyd.warshall.all.pairs.sp}
\title{compute shortest paths for all pairs of nodes}

\description{ compute shortest paths for all pairs of nodes }

\usage{
floyd.warshall.all.pairs.sp(g)
}

\arguments{
  \item{g}{ graph object with edge weights given}
}

\details{
Compute shortest paths between every pair of vertices for a dense graph.  
It works on both undirected and directed graph.
The result is given as a distance matrix.  The matrix is symmetric for an
undirected graph, and asymmetric (very likely) for a directed graph.
For a sparse graph, the \code{johnson.all.pairs.sp} functions
should be used instead.

See documentation on these algorithms in Boost Graph Library for more details.
}

\value{
A matrix of shortest path lengths between all pairs of nodes in the graph.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{johnson.all.pairs.sp}}} 

\examples{
con <- file(system.file("XML/conn.gxl", package="RBGL"), open="r")
coex <- fromGXL(con)
close(con)
floyd.warshall.all.pairs.sp(coex)
}
\keyword{ models }
