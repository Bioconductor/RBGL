\name{brandes.betweenness.centrality}
\alias{brandes.betweenness.centrality}
\title{Compute betweenness centrality for an undirected graph}
\description{Compute betweenness centrality for an undirected graph}
}
\usage{
brandes.betweenness.centrality(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}
\details{
\code{brandes.betweenness.centrality} computes the betweenness centrality of 
each vertex or each edge in the graph.  
The edge betweenness centrality indicates for each edge the betweenness 
centrality that was contributed to the target(s) of the edge (plural for 
undirected graphs).  
Edge betweenness centrality can be used to determine the edges through which 
most shortest paths must pass.
Relative betweenness centrality scales each absolute centrality by 
\dQuote{2/(n*n - 3n + 2)} (where n is the number of vertices in the graph).
Dominance vertex is a measure of the maximum "betweenness" of any point in the graph: it will be 0 for complete graphs and 1 for "wheel" graphs.
}
\value{
\item{betweenness.centrality.vertices}{betweenness centrality of each vertex}
\item{betweenness.centrality.edges}{betweenness centrality of each edge}
\item{relative.betweenness.centrality.vertices}{relative betweenness centrality of each vertex}
\item{dominance}{maximum betweenness of any point in the graph}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex@edgemode <- "undirected"
brandes.betweenness.centrality(coex)
}
\keyword{ models }
