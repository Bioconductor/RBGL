\name{edgeConnectivity}
\alias{edgeConnectivity}
\title{computed edge connectivity and min disconnecting set
for an undirected graph}
\description{computed edge connectivity and min disconnecting set
for an undirected graph}
}
\usage{
edgeConnectivity(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}
\details{
Consider a graph G consisting of a single connected component.
The edge connectivity of G is the minimum number
of edges in G that can be cut to produce a graph with two
(disconnected) components.  The set of edges in this cut
is called the minimum disconnecting set.
}
\value{
A list:
\item{connectivity}{the integer describing the number of edges that must
be severed to obtain two components}
\item{minDisconSet}{a list (of length \code{connectivity}) of pairs of node names describing the edges
that need to be cut to obtain two components}
}
\references{Boost Graph Library by Siek et al.  }
\author{Vince Carey <stvjc@channing.harvard.edu>}

\seealso{\code{\link{minCut}}, \code{\link{edmunds.karp.max.flow}}, \code{\link{push.relabel.max.flow}}  }

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
edgeConnectivity(coex)
}
\keyword{ models }
