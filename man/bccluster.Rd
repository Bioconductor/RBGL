\name{betweenness.centrality.clustering}
\alias{betweenness.centrality.clustering}
\title{Graph clustering based on edge betweenness centrality}
\description{Graph clustering based on edge betweenness centrality}
}
\usage{
betweenness.centrality.clustering(g, threshold = -1, normalize = T))
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
This algorithm implements graph clustering based on edge betweenness centrality.
It is an iterative algorithm, where in each step it compute the edge betweenness
centrality and removes the edge with the maximum betweenness centrality.
}
\value{
\item{edges}{}
\item{betweenness.centrality.clustering}{}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex@edgemode <- "undirected"
betweenness.centrality.clustering(coex)
}
\keyword{ models }
