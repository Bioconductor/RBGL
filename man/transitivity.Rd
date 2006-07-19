\name{transitivity}
\alias{transitivity}
\title{Calculate transitivity for an undirected graph}
\description{Calculate transitivity for an undirected graph }
}
\usage{
transitivity(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class }
}
\details{
For an undirected graph {G}, let delta(v) be the number of triangles with {v} asa node, let tau(v) be the number of paths of length 2 with {v} as the center
node.

Define transitivity T(G) = sum(delta(v)) / sum(tau(v)), for all v in V.
}
\value{
Transitivity for graph \code{g}.
}

\references{
Approximating Clustering Coefficient and Transitivity, T. Schank, D. Wagner,
Journal of Graph Algorithms and Applications, Vol. 9, No. 2 (2005).
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{ClusteringCoef}

\examples{
g <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
tc <- transitivity(g)
}

\keyword{ models }
