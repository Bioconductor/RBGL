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
You can find the definitions of betweenness centrality for vertex and edges and the details of the algorithm in 
the documentation on brandes betweenness centrality in Boost Graph Library.
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
