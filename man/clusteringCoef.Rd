\name{clusteringCoef}
\alias{clusteringCoef}
\title{Calculate clustering coefficient for an undirected graph}
\description{Calculate clustering coefficient for an undirected graph }

\usage{
clusteringCoef(g, Weighted=FALSE, vW=degree(g))
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
  \item{Weighted}{calculate weighted clustering coefficient or not}
  \item{vW}{vertex weights to use when calculating weighted clustering coefficient}
}

\details{
For an undirected graph \code{G}, let delta(v) be the number of triangles with
\code{v}  
as a node, let tau(v) be the number of triples, i.e., paths of length 2 with 
\code{v} as the center node.  

Let V' be the set of nodes with degree at least 2.

Define clustering coefficient for \code{v}, c(v) = (delta(v) / tau(v)). 

Define clustering coefficient for \code{G}, C(G) = sum(c(v)) / |V'|, 
for all \code{v} in V'. 

Define weighted clustering coefficient for \code{g}, 
Cw(G) = sum(w(v) * c(v)) / sum(w(v)), for all \code{v} in V'.
}

\value{
Clustering coefficient for graph \code{G}.
}

\references{
Approximating Clustering Coefficient and Transitivity, T. Schank, D. Wagner, 
Journal of Graph Algorithms and Applications, Vol. 9, No. 2 (2005).
}

\author{Li Long \email{li.long@isb-sib.ch}}

\seealso{clusteringCoefAppr, transitivity, graphGenerator}

\examples{
con <- file(system.file("XML/conn.gxl",package="RBGL"))
g <- fromGXL(con)
close(con)
cc <- clusteringCoef(g)
ccw1 <- clusteringCoef(g, Weighted=TRUE)
vW  <- c(1, 1, 1, 1, 1,1, 1, 1)
ccw2 <- clusteringCoef(g, Weighted=TRUE, vW)
}

\keyword{ models }
