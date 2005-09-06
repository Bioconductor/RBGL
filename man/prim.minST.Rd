\name{prim.minST}
\alias{prim.minST}
\title{Compute minimum spanning tree for an undirected graph}
\description{Compute minimum spanning tree for an undirected graph}
}
\usage{
prim.minST(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}
\details{
This is Prim's algorithm for solving the minimum spanning tree problem for an 
undirected graph with weighted edges. 
See documentations on this function in Boost Graph Library for more details.
}
\value{
  A list of
\item{edges}{the edges that form the minimum spanning tree}
\item{weights}{the total weight of the minimum spanning tree}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{mstree.kruskal}}}

\examples{
coex <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
prim.minST(coex)
}
\keyword{ models }
