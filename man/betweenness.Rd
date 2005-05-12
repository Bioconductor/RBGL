\name{brandesBetweennessCentrality}
\alias{brandesBetweennessCentrality}
\title{Compute betweenness centrality for an undirected graph}
\description{Compute betweenness centrality for an undirected graph}
}
\usage{
brandesBetweennessCentrality(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}
\details{
\code{brandesBetweennessCentrality} computes the betweenness centrality of 
each vertex or each edge in the graph.  
You can find the definitions of betweenness centrality for vertex and edges and the details of the algorithm in 
the documentation on brandes betweenness centrality in Boost Graph Library.
}
\value{
  A list of
\item{betweennessCentrality.vertices}{betweenness centrality of each vertex}
\item{betweennessCentrality.edges}{betweenness centrality of each edge}
\item{relativeBetweennessCentrality.vertices}{relative betweenness centrality of each vertex}
\item{dominance}{maximum betweenness of any point in the graph}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex@edgemode <- "undirected"
brandesBetweennessCentrality(coex)
}
\keyword{ models }
