\name{bfs}
\alias{bfs}
\alias{dfs}
\alias{dfs,graph-method}
\alias{bfs,graph,numeric,logical-method}
\title{ breadth/depth-first search using boost C++ }
\description{
returns integer vector of indices with
order of discovery by BFS or DFS
}
\usage{
bfs(graph, init.node=1, checkConn=FALSE) 
\S4method{dfs}{graph}{graph}
\S4method{bfs}{graph,numeric,logical}{graph,init.node,checkConn}
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
\note{Very preliminary}

\examples{
dd <- fromGXL(file(system.file("XML/bfsex.gxl",package="RBGL")))
dd@edgemode <- "undirected"
bfs(dd)
bfs(dd,2)

dd <- fromGXL(file(system.file("XML/dfsex.gxl",package="RBGL")))
dfs(dd)
}
\keyword{ graphs }
