\name{planarFaceTraversal}
\alias{planarFaceTraversal}

\title{planarFaceTraversal}

\description{planarFaceTraversal description }

\usage{
planarFaceTraversal(g) 
}

\arguments{
  \item{g}{ instance of class graphNEL from Bioconductor graph class}
}

\details{
}

\value{
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

library("RBGL")

V <- LETTERS[1:9]
g <- new("graphNEL", nodes=V, edgemode="undirected")
g <- addEdge(V[1+0],V[1+1],g)
g <- addEdge(V[1+1],V[1+2],g)
g <- addEdge(V[1+3],V[1+4],g)
g <- addEdge(V[1+4],V[1+5],g)
g <- addEdge(V[1+6],V[1+7],g)
g <- addEdge(V[1+7],V[1+8],g)
g <- addEdge(V[1+0],V[1+3],g)
g <- addEdge(V[1+3],V[1+6],g)
g <- addEdge(V[1+1],V[1+4],g)
g <- addEdge(V[1+4],V[1+7],g)
g <- addEdge(V[1+2],V[1+5],g)
g <- addEdge(V[1+5],V[1+8],g)

x1 <- planarFaceTraversal(g)
x1

 }

\keyword{ graphs }

