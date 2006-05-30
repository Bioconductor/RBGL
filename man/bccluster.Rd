\name{betweenness.centrality.clustering}
\alias{betweenness.centrality.clustering}
\title{Graph clustering based on edge betweenness centrality}
\description{Graph clustering based on edge betweenness centrality}
}
\usage{
betweenness.centrality.clustering(g, threshold = -1, normalize = T)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{threshold}{threshold to terminate clustering process}
  \item{normalize}{boolean, when TRUE, the  edge betweenness centrality is 
     scaled by \code{2/((n-1)(n-2))} where \code{n} is the number of vertices 
     in \code{g}; when FALSE, the edge betweenness centrality is the absolute
     value}
}
\details{
To implement graph clustering based on edge betweenness centrality.

The algorithm is iterative, at each step it computes the edge betweenness
centrality and removes the edge with maximum betweenness centrality when it
is above the given \code{threshold}.  When the maximum betweenness centrality 
falls below the threshold, the algorithm terminates.

See documentation on Clustering algorithms in Boost Graph Library for details.
}
\value{
A list of
\item{no.of.edges}{number of remaining edges after removal}
\item{edges}{remaining edges}
\item{edge.betweenness.centrality}{betweenness centrality of remaining edges}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{brandes.betweenness.centrality}}} 

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex <- ugraph(coex)
betweenness.centrality.clustering(coex, 0.1, TRUE)
}
\keyword{ models }
