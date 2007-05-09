\name{mstree.prim}
\alias{mstree.prim}
\alias{prim.minST}

\title{Compute minimum spanning tree for an undirected graph}

\description{Compute minimum spanning tree for an undirected graph}

\usage{
mstree.prim(g)
prim.minST(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}

\details{
This is Prim's algorithm for solving the minimum spanning tree problem for an 
undirected graph with weighted edges. 

The function \code{prim.minST} is depricated and renamed, use 
\code{mstree.prim} instead.

See documentations on this function in Boost Graph Library for more details.
}

\value{
  A list of
\item{edges}{the edges that form the minimum spanning tree}
\item{weights}{the total weight of the minimum spanning tree}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{mstree.kruskal}}}

\examples{
coex <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
mstree.prim(coex)
}
\keyword{ models }
