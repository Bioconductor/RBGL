\name{edmondsOptimumBranching}
\alias{edmondsOptimumBranching}

\title{edmondsOptimumBranching}

\description{edmondsOptimumBranching description }

\usage{
edmondsOptimumBranching(g) 
}

\arguments{
  \item{g}{ instance of class graphNEL from Bioconductor graph class}
}

\details{
This is an implementation of Edmonds' algorithm to find optimum branching in a directed graph.  See references for details.
}

\value{
}

\references{
See Edmonds' Algorithm on http://edmonds-alg.sourceforge.net.
}

\author{ Li Long <li.long@isb-sib.ch> }

\examples{
V <- LETTERS[1:4]
g <- new("graphNEL", nodes=V, edgemode="directed")
g <- addEdge(V[1+0],V[1+1],g, 3)
g <- addEdge(V[1+0],V[2+1],g, 1.5)
g <- addEdge(V[1+0],V[3+1],g, 1.8)
g <- addEdge(V[1+1],V[2+1],g, 4.3)
g <- addEdge(V[1+2],V[3+1],g, 2.2)

x11 <- edmondsOptimumBranching(g)
x11
}

\keyword{ graphs }

