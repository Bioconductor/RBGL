\name{dominatorTree}
\alias{dominatorTree}
\alias{lengauerTarjanDominatorTree}

\title{Compute dominator tree from a vertex in a directed graph }

\description{Compute dominator tree from a vertex in a directed graph}

\usage{
dominatorTree(g, start=nodes(g)[1])
lengauerTarjanDominatorTree(g, start=nodes(g)[1])
}

\arguments{
  \item{g}{a directed graph, one instance of the \code{graph} class }
  \item{start}{a vertex in graph \code{g} }
}

\details{
As stated in documentation on Lengauer Tarjan dominator tree in Boost Graph Library:

A vertex u dominates a vertex v, if every path of directed graph from the entry to v must go through u. 

This function builds the dominator tree for a directed graph.
}

\value{
  Output is a vector, giving each node its immediate dominator.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{ }

\examples{
con1 <- file(system.file("XML/dominator.gxl",package="RBGL"), open="r")
g1 <- fromGXL(con1)
close(con1)

dominatorTree(g1)
lengauerTarjanDominatorTree(g1)
}
\keyword{ models }
