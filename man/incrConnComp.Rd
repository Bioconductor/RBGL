\name{incremental.components}
\alias{incremental.components}
\title{Compute connected components for an undirected graph}
\description{Compute connected components for an undirected graph}
}
\usage{
init.incremental.components(g)
incremental.components(g)
same.component(g, v1, v2)
}
\arguments{
  \item{g}{an instance of the \code{graph} class }
  \item{v1}{one vertex of the given graph }
  \item{v2}{another vertex of the given graph }
}
\details{
This is a family of functions that work together to calculate the connected 
components of an undirected graph. The algorithm used here is based on the 
disjoint-sets (fast union-find) data structure which is a good method to use for
situations where the graph is growing (edges are being added) and the connected 
components information needs to be updated repeatedly. This method does not 
cover the situation where edges are both added and removed from the graph, 
hence it is called incremental (and not fully dynamic). 
Currently, the codes can only handle ONE incremental graph at a time.  When you
start working on another graph by calling "init.incremental.components", the
disjoint-sets info on the previous graph is lost. 
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
