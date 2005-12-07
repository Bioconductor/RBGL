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

The result is given as a distance matrix.

For a sparse matrix, use \code{johnson.all.pairs.sp} instead.

See documentation on these algorithms in Boost Graph Library for more details.
}
\value{
matrix of shortest path lengths, read from row node to col node
}
\references{Boost Graph Library by Siek et al. }
\author{Li Long <li.long@isb-sib.ch>}
\note{  }

\seealso{\code{\link{johnson.all.pairs.sp}}} 

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl", package="RBGL")))
floyd.warshall.all.pairs.sp(coex)
}
\keyword{ models }
