\name{dijkstra.sp}
\alias{dijkstra.sp}
\alias{sp.between}
\alias{sp.between.scalar}
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
  \item{init.ind}{ if numeric, numerical index (1 based) of where to start search;
if character, name of node where search begins}
  \item{start}{ character vector: node name for start of path}
  \item{finish}{ character vector: node name for end of path}
}
\details{
calls to boost dijkstra shortest paths
}
\value{
dijkstra.sp returns list of two vectors of nodes, shortest distance to each node,
and penultimate node of shortest path, and of the index of the
staring node, assumes 1-based counting

If a distance is requested between two nodes that are disconnected,
Inf is returned.

sp.between returns a list whose size depends upon the start and
finish arguments.  If start and finish are both atoms, then
a list is returned with elements the length of the shortest
path and the vector of nodes comprising that path.  If one of start
or finish is an n-vector, n>1, and the other is an atom,
a list is returned with n elements, each element of which is a list
possessing a length and path element.  If both start and finish
are n-vectors, a list of n elements is returned, with ith element the distance
and path between start[i] and finish[i].  It is an error to have
start and finish of different nonunity lengths.
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{Very preliminary}

\examples{
dd <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
dijkstra.sp(dd)
dijkstra.sp(dd,2)
ospf <- fromGXL(file(system.file("XML/ospf.gxl",package="RBGL")))
dijkstra.sp(ospf,6)
dijkstra.sp(ospf,nodes(ospf)[6])
sp.between(ospf, "RT6", "RT1")
sp.between(ospf, c("RT6", "RT2"), "RT1")
sp.between(ospf, c("RT6", "RT2"), c("RT1","RT5"))
if (interactive()) # see an error for query on nonexistent path
sp.between(ospf,"N10", "N13")
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
