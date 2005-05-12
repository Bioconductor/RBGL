\name{dijkstraSP}
\alias{dijkstraSP}
\title{ Dijkstra's shortest paths using boost C++ }
\description{
dijkstra's shortest paths
}
\usage{
dijkstraSP(g,start=nodes(g)[1])
}
\arguments{
  \item{g}{ instance of class graph }
  \item{start}{ character: node name for start of path }
}
\details{
These functions are interfaces to the Boost graph library C++
routines for Dijkstra's shortest paths.
}

\value{
A list with elements:
\item{distance}{The vector of distances from \code{start} to each node
of \code{g}; includes \code{Inf} when there is no path from \code{start}.}
\item{penult}{A vector of indices 
(in \code{nodes(g)}) of predecessors corresponding to each node
on the path from that node back to \code{start}}. For example, if the
element one of this vector has value \code{10}, that means that the
predecessor of node \code{1} is node \code{10}. The next predecessor is
found by examining \code{penult[10]}.
\item{start}{The start node that was supplied in the call to
\code{dijkstraSP}.}
}

\author{ VJ Carey <stvjc@channing.harvard.edu> }

\seealso{\code{\link{spBetween}}}

\examples{
dd <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
dijkstraSP(dd)
dijkstraSP(dd,nodes(dd)[2])
ospf <- fromGXL(file(system.file("XML/ospf.gxl",package="RBGL")))
dijkstraSP(ospf,nodes(ospf)[6])
}
\keyword{ graphs }
