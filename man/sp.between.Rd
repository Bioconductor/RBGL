\name{sp.between}
\alias{sp.between}
\alias{sp.between.old}
\alias{sp.between.scalar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Dijkstra's shortest paths using boost C++ }
\description{
dijkstra's shortest paths
}
\usage{
sp.between.scalar(g,start,finish, eW=unlist(edgeWeights(g))) 
sp.between(g,start,finish) 
sp.between.old(g,start,finish) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{g}{ instance of class graph }
  \item{start}{ character: node name(s) for start of path  }
  \item{finish}{ character vector: node names for end of path}
  \item{eW}{ numeric: edge weights.  }
}
\details{
These functions are interfaces to the Boost graph library C++
routines for Dijkstra's shortest paths.
}

\value{
\code{sp.between.scalar} returns a list with elements:
\item{length}{length (using edge weights) of shortest path}
\item{path}{vector of node names describing the path}
\item{pweights}{named vector of edge weights for all
weights constituting the path -- see \code{\link[graph]{pathWeights}}
for caveats about undirected graph representation.}

Note that if \code{start} or \code{finish} are vectors,
a list of results of repeated application of \code{sp.between.scalar}
is returned.
}


\note{
sp.between returns a list whose size depends upon the start and
finish arguments.  If start and finish are both scalar, then
a list is returned with elements the length of the shortest
path and the vector of nodes comprising that path.  If one of start
or finish is an n-vector, n>1, and the other is an atom,
a list is returned with n elements, each element of which is a list
possessing a length and path element.  If both start and finish
are n-vectors, a list of n elements is returned, with ith element the distance
and path between start[i] and finish[i].  It is an error to have
start and finish of different nonunity lengths.

sp.between.old is included for archival purposes.
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\seealso{\code{\link{bellman.ford.sp}}, \code{\link{dag.sp}}, \code{\link{dijkstra.sp}}, \code{\link{johnson.all.pairs.sp}}}

\examples{
ospf <- fromGXL(file(system.file("XML/ospf.gxl",package="RBGL"), open="r"))
dijkstra.sp(ospf,nodes(ospf)[6])
sp.between(ospf, "RT6", "RT1")
sp.between(ospf, c("RT6", "RT2"), "RT1")
sp.between(ospf, c("RT6", "RT2"), c("RT1","RT5"))
if (interactive()) # see an error for query on nonexistent path
sp.between(ospf,"N10", "N13")
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
