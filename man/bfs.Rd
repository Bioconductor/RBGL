\name{bfs}
\alias{bfs}
\alias{dfs}
\alias{dfs,graph-method}
\alias{bfs,graph-method}
\alias{bfs,graph,ANY,ANY-method}
\title{ Breadth and Depth-first search }
\description{
  These functions carry out breadth and depth first search using
  routines from the BOOST library. For \code{dfs} a \code{list}
  containing two \code{vectors} is returned. These indicate the sequence
  in which a node is first visited and last visited. For \code{dfs} a
  single vector is returned indicating the order in which a node was visited.
}
\usage{
bfs(graph, init.node=1, checkConn=FALSE)
\S4method{dfs}{graph}(graph)
\S4method{bfs}{graph}(graph)
}
\arguments{
  \item{graph}{ instance of class graph from Bioconductor graph class}
  \item{init.node}{ index (1 based) of where to start search}
  \item{checkConn}{ logical indicating whether connectivity of
input graph should be checked }
}
\details{
calls to boost bfs or dfs procedures
}
\value{
for bfs: vector of indices in order of BFS visit

for dfs:
list of two vectors of nodes, order of DFS discovery,
and order of DFS completion (finish)
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\examples{
dd <- fromGXL(file(system.file("XML/bfsex.gxl",package="RBGL")))
dd@edgemode <- "undirected"
bfs(dd)
bfs(dd,2)

dd <- fromGXL(file(system.file("XML/dfsex.gxl",package="RBGL")))
dfs(dd)
}
\keyword{ graphs }
