\name{betweennessCentralityClustering}
\alias{betweennessCentralityClustering}
\title{Graph clustering based on edge betweenness centrality}
\description{Graph clustering based on edge betweenness centrality}
}
\usage{
betweennessCentralityClustering(g, threshold = -1, normalize = T))
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{threshold}{threshold to terminate clustering process}
  \item{normalize}{boolean, when true, the threshold is compared with the 
     normalized edge centrality based on the input graph; when false, the 
     threshold is compared with the absolute edge centrality}
}
\details{
To implement graph clustering based on edge betweenness centrality.
The algorithm is iterative, at each step it computes the edge betweenness
centrality and removes the edge with the maximum betweenness centrality.
See documentation on Clustering algorithms in Boost Graph Library for details.
}
\value{
A list of
\item{betweennessCentralityVertices}{betweenness centrality of each vertex}
\item{betweennessCentralityEdges}{betweenness centrality of each edge}
\item{relativeBetweennessCentralityVertices}{relative betweenness centrality of each vertex}
\item{dominance}{central point dominance}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex@edgemode <- "undirected"
betweennessCentralityClustering(coex)
}
\keyword{ models }
