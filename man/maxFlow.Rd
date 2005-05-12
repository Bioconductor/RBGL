\name{maxFlow}
\alias{edmundsKarpMaxFlow}
\alias{pushRelabelMaxFlow}
\title{Compute max flow for a directed graph }
\description{Compute max flow for a directed graph }
}
\usage{
edmundsKarpMaxFlow(g, source, sink)
pushRelabelMaxFlow(g, source, sink)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{directed}}
  \item{source}{ node name (character) or node number (int) for the source of the flow }
  \item{sink}{ node name (character) or node number (int) for the sink of the flow }
}
\details{
Given a directed graph G=(V, E) of a single connected component with a vertex 
\code{source} and a vertex \code{sink}.  Each arc has a positive real valued 
capacity, currently it's equivalent to the weight of the arc.  The flow of the 
network is the net flow entering the vertex \code{sink}. The maximum flow 
problem is to determine the maximum possible value for the flow to the 
\code{sink} and the corresponding flow values for each arc. 
See documentation on these algorithms in Boost Graph Library for more details.
}
\value{
  A list of
\item{maxflow}{the max flow from \code{source} to \code{sink} }
\item{edges}{the nodes of the arcs with non-zero capacities}
\item{flows}{the flow values of the arcs with non-zero capacities}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{ }

\examples{
g <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))

ans1 <- edmundsKarpMaxFlow(g, "B", "D")
ans2 <- edmundsKarpMaxFlow(g, 3, 2)     # 3 and 2 equivalent to "C" and "B"

ans3 <- pushRelabelMaxFlow(g, 2, 4)     # 2 and 4 equivalent to "B" and "D"
ans4 <- pushRelabelMaxFlow(g, "C", "B")
}
\keyword{ models }
