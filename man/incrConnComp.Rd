\name{incremental.components}
\alias{incremental.components}
\title{Compute connected components for an undirected graph}
\description{Compute connected components for an undirected graph}
}
\usage{
init.incremental.components(g)
incremental.components(g)
same.component(g, node1, node2)
}
\arguments{
  \item{g}{an instance of the \code{graph} class }
  \item{node1}{one vertex of the given graph }
  \item{node2}{another vertex of the given graph }
}
\details{
This family of functions work together to calculate the connected 
components of an undirected graph. The algorithm is based on the 
disjoint-sets.  It works where the graph is growing by adding new edges.
Call "init.incremental.components" to initialize the calculation on a new graph.
Call "incremental.components" to re-calculate connected components after 
growing the graph.  Call "same.component" to learn if two given vertices 
are in the same connected components. 
Currently, the codes can only handle ONE incremental graph at a time.  When you
start working on another graph by calling "init.incremental.components", the
disjoint-sets info on the previous graph is lost. 
See documentation on Incremental Connected Components in Boost Graph Library 
for more details.
}
\value{
\item{}{the number of connected components in the current graph }
\item{}{connected components in the current graph }
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{  }

\examples{
coex <- fromGXL(file(system.file("XML/conn2.gxl",package="RBGL")))
init.incremental.components(coex)
incremental.components(coex)
v1 <- 1
v2 <- 5
same.component(coex, v1, v2)
}
\keyword{ models }
