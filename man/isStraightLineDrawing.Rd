\name{isStraightLineDrawing}
\alias{isStraightLineDrawing}

\title{isStraightLineDrawing}

\description{isStraightLineDrawing description }

\usage{
isStraightLineDrawing(g, drawing) 
}

\arguments{
  \item{g}{ instance of class graphNEL from Bioconductor graph class}
  \item{drawing}{ coordinates of node positions }
}

%\details{
%}
%
%\value{
%}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ Li Long <li.long@isb-sib.ch> }

\examples{

V <- LETTERS[1:7]
g <- new("graphNEL", nodes=V, edgemode="undirected")
g <- addEdge(V[1+0], V[1+1], g)
g <- addEdge(V[1+1], V[2+1], g)
g <- addEdge(V[1+2], V[3+1], g)
g <- addEdge(V[1+3], V[0+1], g)
g <- addEdge(V[1+3], V[4+1], g)
g <- addEdge(V[1+4], V[5+1], g)
g <- addEdge(V[1+5], V[6+1], g)
g <- addEdge(V[1+6], V[3+1], g)
g <- addEdge(V[1+0], V[4+1], g)
g <- addEdge(V[1+1], V[3+1], g)
g <- addEdge(V[1+3], V[5+1], g)
g <- addEdge(V[1+2], V[6+1], g)
g <- addEdge(V[1+1], V[4+1], g)
g <- addEdge(V[1+1], V[5+1], g)
g <- addEdge(V[1+1], V[6+1], g)


x3 <- chrobakPayneStraightLineDrawing(g)

x8 <- isStraightLineDrawing(g, x3)
x8

 }

\keyword{ graphs }

