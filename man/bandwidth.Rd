\name{bandwidth}
\alias{bandwidth}
\title{Compute bandwidth for an undirected graph}
\description{Compute bandwidth for an undirected graph}
}
\usage{
bandwidth(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}
\details{
The bandwidth of an undirected graph G=(V, E) is the maximum distance between 
two adjacent vertices, with distance measured on a line upon which the vertices 
have been placed at unit intervals. To put it another way, if the vertices of 
an undirected graph G=(V,E) are each assigned an index from zero to |V| - 1 
given by index[v], then the bandwidth of G is

B(G) = max { |index[u] - index[v]|  | (u,v) in E }
}
\value{
\item{bandwidth}{the bandwidth of the given graph}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{  }

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
coex@edgemode <- "undirected"
bandwidth(coex)
}
\keyword{ models }
