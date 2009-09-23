\name{extractPath}
\alias{extractPath}
\title{convert a dijkstra.sp predecessor structure into the path joining two nodes   }

\description{
determine a path between two nodes in a graph,
using output of \code{\link{dijkstra.sp}}.
}

\usage{
extractPath(s, f, pens)
}

\arguments{
  \item{s}{ index of starting node in nodes vector
of the graph from which \code{pens} was derived}
  \item{f}{ index of ending node in nodes vector }
  \item{pens}{ predecessor index vector as returned
in the \code{preds} component of \code{\link{dijkstra.sp}} output}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Vince Carey <stvjc@channing.harvard.edu>}

\seealso{  \code{\link[e1071:shortestPaths]{allShortestPaths}}}

\examples{
data(FileDep)
dd <- dijkstra.sp(FileDep)
extractPath(1,9,dd$pen)
}
\keyword{ models }
