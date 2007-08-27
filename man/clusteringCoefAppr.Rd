\name{clusteringCoefAppr}
\alias{clusteringCoefAppr}
\title{Approximate clustering coefficient for an undirected graph}
\description{Approximate clustering coefficient for an undirected graph }

\usage{
clusteringCoefAppr(g, k=length(nodes(g)), Weighted=FALSE, vW=degree(g))
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
  \item{Weighted}{calculate weighted clustering coefficient or not}
  \item{vW}{vertex weights to use when calculating weighted clustering coefficient}
  \item{k}{parameter controls total expected runtime}
}

\details{
It is quite expensive to compute cluster coefficient and transitivity exactly 
for a large graph by computing the number of triangles in the graph.  Instead, 
\code{clusteringCoefAppr} samples triples with appropriate probability, returns
the ratio between the number of existing edges and the number of samples.   

MORE ABOUT CHOICE OF K.

See reference for more details.
}

\value{
Approximated clustering coefficient for graph \code{g}.
}

\references{
Approximating Clustering Coefficient and Transitivity, T. Schank, D. Wagner, 
Journal of Graph Algorithms and Applications, Vol. 9, No. 2 (2005).
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{clusteringCoef, transitivity, graphGenerator}

\examples{
con <- file(system.file("XML/conn.gxl",package="RBGL"))
g <- fromGXL(con)
close(con)

k = length(nodes(g))
cc <- clusteringCoefAppr(g, k)
ccw1 <- clusteringCoefAppr(g, k, Weighted=TRUE)
vW  <- c(1, 1, 1, 1, 1,1, 1, 1)
ccw2 <- clusteringCoefAppr(g, k, Weighted=TRUE, vW)
}

\keyword{ models }
