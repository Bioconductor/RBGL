\name{dijkstra.sp}
\alias{dijkstra.sp}
%\alias{sp.between}
%\alias{sp.between.scalar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Dijkstra's shortest paths using boost C++ }
\description{
dijkstra's shortest paths
}
\usage{
dijkstra.sp(g,start=nodes(g)[1])
%sp.between(g,start,finish) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{g}{ instance of class graph }
  \item{start}{ character: node name for start of path }
%must have length 1 for dijkstra.sp}
%  \item{finish}{ character vector: node name for end of path}
}
\details{
These functions are interfaces to the Boost graph library C++
routines for Dijkstra's shortest paths.
}

\value{
a list with elements:
\item{distance}{vector of distances from \code{start} to each node
of \code{g}; includes \code{Inf} when there is no path from \code{start}}
\item{pred}{vector of indices 
(in \code{nodes(g)}) of predecessors corresponding to each node
on the path from that node back to \code{start}}
}

\author{ VJ Carey <stvjc@channing.harvard.edu> }

\seealso{\code{\link{sp.between}}}

\examples{
dd <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
dijkstra.sp(dd)
dijkstra.sp(dd,nodes(dd)[2])
ospf <- fromGXL(file(system.file("XML/ospf.gxl",package="RBGL")))
dijkstra.sp(ospf,nodes(ospf)[6])
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
