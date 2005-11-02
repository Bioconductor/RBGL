\name{minCut}
\alias{minCut}
\title{Compute min-cut for an undirected graph}
\description{Compute min-cut for an undirected graph}
}
\usage{
minCut(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}
\details{
Given an undirected graph G=(V, E) of a single connected component, a cut is 
a partition of the set of vertices into two non-empty subsets S and V-S, a 
cost is the number of edges that are incident on one vertex in S and one 
vertex in V-S.  The min-cut problem is to find a cut (S, V-S) of minimum cost.

For simplicity, the returned subset S is the smaller of the two subsets.

}
\value{
  A list of 
\item{mincut}{the number of edges to be severed to obtain the minimum cut}
\item{S}{the smaller subset of vertices in the minimum cut}
\item{V-S}{the other subset of vertices in the minimum cut}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{edgeConnectivity}}}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
minCut(coex)
}
\keyword{ models }
