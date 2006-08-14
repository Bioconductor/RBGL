\name{graphGenerator}
\alias{graphGenerator}
\title{Generate an undirected graph with adjustable clustering coefficient}
\description{Generate an undirected graph with adjustable clustering coefficient }
}
\usage{
graphGenerator(n, d, o)
}
\arguments{
  \item{n}{no. of nodes in the generated graph}
  \item{d}{parameter for preferential attachment}
  \item{o}{parameter for triple generation}
}
\details{
The graph generator works according to the prefential attachment rule.  It also generates graphs with adjustable clustering coefficient.  Parameter \code{d} specifies how many preferred edges a new node has.  Parameter \code{o} limits how many triples to add to a new node.

See reference for details.
}
\value{
\item{no. of nodes}{No. of nodes in the generated graph}
\item{no. of edges}{No. of edges in the generated graph}
\item{edges}{Edges in the generated graph}
}

\references{
Approximating Clustering Coefficient and Transitivity, T. Schank, D. Wagner, 
Journal of Graph Algorithms and Applications, Vol. 9, No. 2 (2005).
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{clusteringCoef, transitivity, clusteringCoefAppr}

\examples{
n <- 20
d <- 6
o <- 3
gg <- graphGenerator(n, d, o)
}

\keyword{ models }
