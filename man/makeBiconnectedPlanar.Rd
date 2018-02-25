\name{makeBiconnectedPlanar}
\alias{makeBiconnectedPlanar}

\title{makeBiconnectedPlanar}

\description{makeBiconnectedPlanar description }

\usage{
makeBiconnectedPlanar(g) 
}

\arguments{
  \item{g}{ instance of class graphNEL from Bioconductor graph class}
}

\details{
From 

http://www.boost.org/doc/libs/1_49_0/libs/graph/doc/planar_graphs.html


An undirected graph is connected if, for any two vertices 
u and v, there's a path from u to v. 
An undirected graph is biconnected if it is connected and 
it remains connected even if any single vertex is removed. 
Finally, a planar graph is maximal planar (also called triangulated) 
if no additional edge (with the exception of self-loops and 
parallel edges) can be added to it without creating a 
non-planar graph. Any maximal planar simple graph on n > 
2 vertices has exactly 3n - 6 edges and 2n 
- 4 faces, a consequence of Euler's formula. If a 
planar graph isn't connected, isn't biconnected, or isn't maximal planar, 
there is some set of edges that can be added 
to the graph to make it satisfy any of those 
three properties while preserving planarity. Many planar graph drawing algorithms 
make at least one of these three assumptions about the 
input graph, so there are functions in the Boost Graph 
Library that can help:

make_connected adds a minimal set 
of edges to an undirected graph to make it connected.

make_biconnected_planar adds a set of edges to a connected, undirected 
planar graph to make it biconnected while preserving planarity.

make_maximal_planar 
adds a set of edges to a biconnected, undirected planar 
graph to make it maximal planar.

The function documented here implements the second approach.
}

\value{
A list with two elements: `Is planar` is a logical indicating
achievement of planarity, and `new graph`, a graphNEL instance
that is biconnected and planar.
}


\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ Li Long <li.long@isb-sib.ch> }

\examples{

V <- LETTERS[1:11]
g <- new("graphNEL", nodes=V, edgemode="undirected")
g <- addEdge(V[1+0], V[1+1], g)
g <- addEdge(V[1+2], V[3+1], g)
g <- addEdge(V[1+3], V[0+1], g)
g <- addEdge(V[1+3], V[4+1], g)
g <- addEdge(V[1+4], V[5+1], g)
g <- addEdge(V[1+5], V[3+1], g)
g <- addEdge(V[1+5], V[6+1], g)
g <- addEdge(V[1+6], V[7+1], g)
g <- addEdge(V[1+7], V[8+1], g)
g <- addEdge(V[1+8], V[5+1], g)
g <- addEdge(V[1+8], V[9+1], g)
g <- addEdge(V[1+0], V[10+1], g)

x6 <- makeBiconnectedPlanar(g)
x6

 }

\keyword{ graphs }

