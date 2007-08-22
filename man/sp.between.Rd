\name{sp.between}
\alias{sp.between}
\alias{sp.between.scalar}

\title{ Dijkstra's shortest paths using boost C++ }

\description{ dijkstra's shortest paths }

\usage{
sp.between(g,start,finish, detail=TRUE) 
}

\arguments{
  \item{g}{ instance of class graph }
  \item{start}{ node name(s) for start of path(s)  }
  \item{finish}{ node name(s) for end of path(s) }
  \item{detail}{ if TRUE, output additional info on the shortest path }
}

\details{
These functions are interfaces to the Boost graph library C++
routines for Dijkstra's shortest paths.

Function \code{sp.between.scalar} is obsolete.  
}

\value{
When \code{start} and/or \code{finish} are vectors, we use the normal cycling 
rule in R to match both vectors and try to find the shortest path for each 
pair.  

Function \code{sp.between} returns a list of info on the shortest paths.  Each
such shortest path is designated by its starting node and its ending node.  
Each element in the returned list contains:
\item{length}{total length (using edge weights) of this shortest path},
\item{path_detail}{if requested, names of the nodes on the shortest path},
\item{length_detail}{if requested, edge weights of this shortest path}.

See \code{\link[graph]{pathWeights}} for caveats about undirected graph representation.
}

\author{ VJ Carey <stvjc@channing.harvard.edu>, Li Long <li.long@isb-sib.ch> }

\seealso{\code{\link{bellman.ford.sp}}, \code{\link{dag.sp}}, \code{\link{dijkstra.sp}}, \code{\link{johnson.all.pairs.sp}}}

\examples{

ospf <- fromGXL(file(system.file("XML/ospf.gxl",package="RBGL"), open="r"))
dijkstra.sp(ospf,nodes(ospf)[6])

sp.between(ospf, "RT6", "RT1")

sp.between(ospf, c("RT6", "RT2"), "RT1", detail=FALSE)

sp.between(ospf, c("RT6", "RT2"), c("RT1","RT5"))

# see NAs for query on nonexistent path
sp.between(ospf,"N10", "N13")

}

\keyword{ graphs }% at least one, from doc/KEYWORDS

