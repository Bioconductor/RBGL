\name{dijkstra.sp}
\alias{dijkstra.sp}
\alias{sp.between}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Dijkstra's shortest paths using boost C++ }
\description{
dijkstra's shortest paths
}
\usage{
dijkstra.sp(x,init.ind=1) # now x assumed to be Bioconductor graph graphNEL
sp.between(g,start,finish) # now x assumed to be Bioconductor graph graphNEL
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ instance of class graphNEL from Bioconductor graph class}
  \item{g}{ instance of class graphNEL from Bioconductor graph class}
  \item{init.ind}{ index (1 based) of where to start search}
  \item{start}{ character atom: node name for start of path}
  \item{finish}{ character atom: node name for end of path}
}
\details{
calls to boost dijkstra shortest paths
}
\value{
list of two vectors of nodes, shortest distance to each node,
and penultimate node of shortest path, and of the
staring node, assumes
1-based counting
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{Very preliminary}

\examples{
dd <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
dijkstra.sp(dd)
dijkstra.sp(dd,2)
ospf <- fromGXL(file(system.file("XML/ospf.gxl",package="RBGL")))
dijkstra.sp(ospf,6)
sp.between(ospf, "RT6", "RT1")
if (interactive()) # see an error for query on nonexistent path
sp.between(ospf,"N10", "N13")
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
