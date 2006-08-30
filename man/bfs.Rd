\name{bfs}
\alias{bfs}
\alias{dfs}

\alias{dfs,graph,character,ANY-method}

\alias{bfs,graph,character,logical-method}
\alias{dfs,graph,character,logical-method}

\alias{bfs,graph,missing,logical-method}
\alias{bfs,graph,character,missing-method}
\alias{dfs,graph,character,missing-method}

\alias{bfs,graph,missing,missing-method}
\alias{dfs,graph,missing,missing-method}

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
bfs(object, node, checkConn=TRUE)
dfs(object, node, checkConn=TRUE)
}
\arguments{
  \item{object}{ instance of class graph from Bioconductor graph class}
  \item{node}{ node name where search starts; defaults to the
node in first position in the node vector.}
  \item{checkConn}{ logical for backwards compatibility; this
parameter has no effect as of RBGL 1.7.9 and will be removed in future versions.}
}
\details{
  These two functions are interfaces to the
  BOOST graph library functions for breadth first
and depth first search. Both methods handle
  unconnected graphs by applying the algorithms
over the connected components. 

  Cormen et al note (p 542) that `results of depth-first search
  may depend upon the order in which the vertices are examined ...
  These different visitation orders tend not to cause problems
  in practice, as any DFS result can usually be used effectively,
  with essentially equivalent results'.  
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
dd <- fromGXL(file(system.file("XML/bfsex.gxl",package="RBGL"), open="r"))
bfs(dd, "r", FALSE)
bfs(dd, "s", TRUE)

dd2 <- fromGXL(file(system.file("XML/dfsex.gxl",package="RBGL"), open="r"))
dfs(dd2, "u", FALSE)
}
\keyword{ graphs }
