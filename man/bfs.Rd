\name{bfs}
\alias{bfs}
\alias{dfs2}
\alias{dfs2,graph-method}
\title{ breadth/depth-first search using boost C++ }
\description{
returns integer vector of indices with
order of discovery by BFS or DFS
}
\usage{
bfs(x, init.ind=1) 
\S4method{dfs2}{graph}{object}
}
\arguments{
  \item{x}{ instance of class graph from Bioconductor graph class}
  \item{object}{ instance of class graph from Bioconductor graph class}
  \item{init.ind}{ index (1 based) of where to start search}
}
\details{
calls to boost bfs
}
\value{
for bfs: vector of indices in order of BFS visit

for dfs2:
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
dfs2(dd)
}
\keyword{ graphs }
