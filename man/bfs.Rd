\name{bfs}
\alias{bfs}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ depth-first search using boost C++ }
\description{
returns integer vector of indices with
order of discovery by BFS
}
\usage{
bfs(x) # now x assumed to be Bioconductor graph graphNEL
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ instance of class graphNEL from Bioconductor graph class}
  \item{init.ind}{ index (1 based) of where to start search}
}
\details{
calls to boost bfs
}
\value{
list of two vectors of nodes, see above, assumes
1-based counting
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{Very preliminary}

\examples{
dd <- fromGXL(file(system.file("XML/bfsex.gxl",package="RBGL")))
dd@edgemode <- "undirected"
bfs(dd)
bfs(dd,2)
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
