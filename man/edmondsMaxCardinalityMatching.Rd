\name{edmondsMaxCardinalityMatching}
\alias{edmondsMaxCardinalityMatching}

\title{edmondsMaxCardinalityMatching}

\description{edmondsMaxCardinalityMatching description }

\usage{
edmondsMaxCardinalityMatching(g) 
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
V <- LETTERS[1:18]
g <- new("graphNEL", nodes=V, edgemode="undirected")
g <- addEdge(V[1+0], V[4+1], g);
g <- addEdge(V[1+1], V[5+1], g);
g <- addEdge(V[1+2], V[6+1], g);
g <- addEdge(V[1+3], V[7+1], g);
g <- addEdge(V[1+4], V[5+1], g);
g <- addEdge(V[1+6], V[7+1], g);
g <- addEdge(V[1+4], V[8+1], g);
g <- addEdge(V[1+5], V[9+1], g);
g <- addEdge(V[1+6], V[10+1], g);
g <- addEdge(V[1+7], V[11+1], g);
g <- addEdge(V[1+8], V[9+1], g);
g <- addEdge(V[1+10], V[11+1], g);
g <- addEdge(V[1+8], V[13+1], g);
g <- addEdge(V[1+9], V[14+1], g);
g <- addEdge(V[1+10], V[15+1], g);
g <- addEdge(V[1+11], V[16+1], g);
g <- addEdge(V[1+14], V[15+1], g);

x9 <- edmondsMaxCardinalityMatching(g)
x9

g <- addEdge(V[1+12], V[13+1], g);
g <- addEdge(V[1+16], V[17+1], g);

x10 <- edmondsMaxCardinalityMatching(g)
x10

 }

\keyword{ graphs }

