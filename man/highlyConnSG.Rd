\name{highlyConnSG}
\alias{highlyConnSG}
\title{Compute highly connected subgraphs for an undirected graph}
\description{Compute highly connected subgraphs for an undirected graph}
}
\usage{
highlyConnSG(g, sat, ldv)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{sat}{singleton adoption threshold, positive integer }
  \item{ldv}{heuristics to remove lower degree vertice, a decreasing sequence of positive integer }
}
\details{
A graph G with n vertices is highly connected if its connectivity k(G) > n/2.  The HCS algorithm partitions a given graph into a set of highly connected subgraphs, by using minimum-cut algorithm recursively.  To improve performance, the approach is refined by adopting singletons, removing low degree vertices and merging clusters.    
}
\value{
  A list of clusters, each is given as vertices in the graph.
}
\references{ A Clustering Algorithm based on Graph Connectivity by E. Hartuv, R. Shamir, 1999.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{edgeConnectivity}}, \code{\link{minCut}}  }

\examples{
coex <- fromGXL(file(system.file("XML/hcs.gxl",package="RBGL")))
highlyConnSG(coex)
}
\keyword{ models }
