\name{bfs}
\alias{bfs}
\alias{dfs}
\alias{dfs,graph,character,ANY-method}
\alias{bfs,graph,character,ANY-method}
\alias{dfs,graph,missing,missing-method}
\alias{bfs,graph,missing,missing-method}
\alias{bfs,graph-method}
\alias{bfs,graph,character-method}
\alias{dfs,graph,character-method}
\alias{bfs,graph,ANY,ANY-method}
\title{ Breadth and Depth-first search }
\description{
  These functions return information on
  graph traversal by breadth and depth first search using
  routines from the BOOST library.  
}
\usage{
bfs(object, node, checkConn=FALSE)
\S4method{dfs}{object,node,checkConn}(graph)
\S4method{bfs}{object,node,checkConn}(graph)
}
\arguments{
  \item{object}{ instance of class graph from Bioconductor graph class}
  \item{node}{ node name where search starts -- note that this is ignored
for dfs }
  \item{checkConn}{ logical indicating whether connectivity of
input graph should be checked }
}
\details{
  These two functions are interfaces to the
  BOOST graph library functions for breadth first
and depth first search. Both methods can handle
  unconnected graphs and unless \code{checkConn} is specified neither
  tests connectivity. If \code{checkConn} is specified and the
  graph is not connected an error is signalled.

  Note that the parameter sequences for bfs and dfs agree, but that
  at present the \code{node} argument to \code{dfs} is ignored.
  Cormen et al note (p 542) that `results of depth-first search
  may depend upon the order in which the vertices are examined ...
  These different visitation orders tend not to cause problems
  in practices, as any DFS result can usually be used effectively,
  with essentially equivalent results'.  To fix the
  starting node of \code{dfs}, order the graph vertex list so that the
  desired start point is the first entry.
}
\value{
For \code{bfs} a vector of node indices in order of BFS visit.

For \code{dfs} a
list of two vectors of nodes, with elements
\code{discover} (order of DFS discovery),
and \code{finish}
(order of DFS completion).
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\examples{
dd <- fromGXL(file(system.file("XML/bfsex.gxl",package="RBGL")))
bfs(dd, "r", FALSE)
bfs(dd, "s", TRUE)

dd2 <- fromGXL(file(system.file("XML/dfsex.gxl",package="RBGL")))
dfs(dd2, "u", FALSE)
}
\keyword{ graphs }
