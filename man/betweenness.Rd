\name{brandes.betweenness.centrality}
\alias{brandes.betweenness.centrality}
\title{Compute betweenness centrality for an undirected graph}
\description{Compute betweenness centrality for an undirected graph}

\usage{
brandes.betweenness.centrality(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}

\details{
\code{Brandes.betweenness.centrality} computes the betweenness centrality of 
each vertex or each edge in the graph, using an algorithm by U. Brandes.

Betweenness centrality of a vertex \code{v} is calculated as follows:
\code{N_st(v)} = no. of shortest paths from \code{s} to \code{t} that pass through \code{v},
\code{N_st} = no. of shortest paths from \code{s} to \code{t},
betweenness centrality of \code{v} = \code{sum(N_st(v)/N_st)} for all vertices \code{s} != \code{v} != \code{t}.

Betweenness centrality of an edge is calculated similarly.

The relative betweenness centrality for a vertex is to scale the betweenness
centrality of the given vertex by \code{2/(n**2 - 3n + 2)} where \code{n} is 
the no. of vertices in the graph.

Central point dominance measures the maximum betweenness of any vertex 
in the graph.

See documentation on brandes betweenness centrality in Boost Graph Library 
for more details.
}

\value{
  A list of
\item{betweenness.centrality.vertices}{betweenness centrality of each vertex}
\item{betweenness.centrality.edges}{betweenness centrality of each edge}
\item{relative.betweenness.centrality.vertices}{relative betweenness centrality of each vertex}
\item{dominance}{maximum betweenness of any point in the graph}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{betweenness.centrality.clustering}}}

% \examples{
% con <- file(system.file("XML/conn.gxl",package="RBGL"), open="r")
% coex <- fromGXL(con)
% close(con)
% coex <- ugraph(coex)
% brandes.betweenness.centrality(coex)
% }
\keyword{ models }

