\name{transitivity}
\alias{transitivity}

\title{Calculate transitivity for an undirected graph}

\description{Calculate transitivity for an undirected graph }

\usage{
transitivity(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
}

\details{
For an undirected graph \code{G}, let delta(v) be the number of triangles with
\code{v}  
as a node, let tau(v) be the number of triples, i.e., paths of length 2 with 
\code{v} as the center node.

Define transitivity T(G) = sum(delta(v)) / sum(tau(v)), for all v in V.
}

\value{
Transitivity for graph \code{g}.
}

\references{
Approximating Clustering Coefficient and Transitivity, T. Schank, D. Wagner,
Journal of Graph Algorithms and Applications, Vol. 9, No. 2 (2005).
}

\author{Li Long \email{li.long@isb-sib.ch}}

\seealso{clusteringCoef, clusteringCoefAppr, graphGenerator}

\examples{
con <- file(system.file("XML/conn.gxl",package="RBGL"))
g <- fromGXL(con)
close(con)

tc <- transitivity(g)
}

\keyword{ models }
