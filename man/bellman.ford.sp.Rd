\name{bellman.ford.sp}
\alias{bellman.ford.sp}
\title{ Bellman-Ford shortest paths using boost C++ }

\description{
Algorithm for the single-source shortest paths problem for a graph with both positive and negative edge weights. 
}

\usage{
bellman.ford.sp(g,start=nodes(g)[1])
}

\arguments{
  \item{g}{ instance of class graph }
  \item{start}{ character: node name for start of path }
}

\details{
This function interfaces to the Boost graph library C++
routines for Bellman-Ford shortest paths.  
Choose the appropriate algorithm to calculate the shortest path 
carefully based on the properties of the given graph.
See documentation on Bellman-Ford algorithm in Boost Graph Library 
for more details.
}

\value{
A list with elements:
\item{all edges minimized}{true if all edges are minimized, false otherwise. }
\item{distance}{The vector of distances from \code{start} to each node
of \code{g}; includes \code{Inf} when there is no path from \code{start}.}
\item{penult}{A vector of indices 
(in \code{nodes(g)}) of predecessors corresponding to each node
on the path from that node back to \code{start}}. For example, if the
element one of this vector has value \code{10}, that means that the
predecessor of node \code{1} is node \code{10}. The next predecessor is
found by examining \code{penult[10]}.
\item{start}{The start node that was supplied in the call to
\code{bellman.ford.sp}.}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ Li Long <li.long@isb-sib.ch> }

\seealso{\code{\link{dag.sp}}, \code{\link{dijkstra.sp}}, \code{\link{johnson.all.pairs.sp}}, \code{\link{sp.between}} }

\examples{
dd <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL"), open="r"))
bellman.ford.sp(dd)
bellman.ford.sp(dd,nodes(dd)[2])
}
\keyword{ graphs }
