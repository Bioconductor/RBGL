\name{kCores}
\alias{kCores}
\title{Find all the k-cores in a graph}

\description{Find all the k-cores in a graph }

\usage{
kCores(g, EdgeType=c("in", "out"))
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
  \item{EdgeType} {what types of edges to be considered when \code{g} is directed }
}

\details{
A k-core in a graph is a subgraph where each node is adjacent to at least a minimum number, k, of the other nodes in the subgraph.

The implementation is based on the algorithm by V. Batagelj and M. Zaversnik, 2002.

The example \code{snacoreex.gxl} is in the paper by V. Batagelj and M. Zaversnik, 2002.
}
\value{
A vector of the core numbers for all the nodes in \code{g}.
}

\references{
Social Network Analysis: Methods and Applications.  By S. Wasserman and K. Faust, pp. 266. 
An O(m) Algorithm for Cores decomposition of networks, by V. Batagelj and M. Zaversnik, 2002.
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
kcoex <- fromGXL(file(system.file("XML/snacoreex.gxl",package="RBGL")))
kCores(kcoex)

kcoex2 <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
kCores(kcoex2)
kCores(kcoex2, "in")
kCores(kcoex2, "out")

}

\keyword{ models }
