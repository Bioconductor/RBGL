\name{dag.sp}
\alias{dag.sp}
\title{ DAG shortest paths using boost C++ }
\description{
Algorithm for the single-source shortest-paths problem on a weighted, directed acyclic graph (DAG)
}
\usage{
dag.sp(g,start=nodes(g)[1])
}
\arguments{
  \item{g}{ instance of class graph }
  \item{start}{ source node for start of paths }
}
\details{
These functions are interfaces to the Boost graph library C++ routines
for single-source shortest-paths on a weighted directed acyclic graph.
Choose appropriate shortest-path algorithms carefully based on the properties
of the input graph.  See documentation in Boost Graph Library for more details.
}

\value{
A list with elements:
\item{distance}{The vector of distances from \code{start} to each node
of \code{g}; includes \code{Inf} when there is no path from \code{start}.}

\item{penult}{A vector of indices
(in \code{nodes(g)}) of predecessors corresponding to each node
on the path from that node back to \code{start}. For example, if the
element one of this vector has value \code{10}, that means that the
predecessor of node \code{1} is node \code{10}. The next predecessor is
found by examining \code{penult[10]}.}

\item{start}{The start node that was supplied in the call to
\code{dag.sp}.}
}

\author{ Li Long <li.long@isb-sib.ch> }

\seealso{}

\examples{
dd <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
dag.sp(dd)
dag.sp(dd,nodes(dd)[2])
}
\keyword{ graphs }

