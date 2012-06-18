\name{biConnComp}
\alias{biConnComp}
\alias{articulationPoints}
\title{Compute biconnected components for a graph}
\description{Compute biconnected components for a graph}

\usage{
biConnComp(g)
articulationPoints(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class} 
}

\details{

A biconnected graph is a connected graph that remains connected when 
any one of its vertices, and all the edges incident on this vertex,
is removed and the graph remains connected.
A biconnected component of a graph is a subgraph which is biconnected.  An 
integer label is assigned to each edge to indicate which biconnected component 
it's in.

A vertex in a graph is called an articulation point if removing it increases 
the number of connected components.

See the documentation for the Boost Graph Library for more details.
}

\value{
For \code{biConnComp}:
a vector whose length is no. of biconnected components, each entry is a list 
of nodes that are on the same biconnected components.

For \code{articulationPoints}:
a vector of articulation points in the graph.

}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

%\seealso{}

\examples{
con <- file(system.file("XML/conn.gxl",package="RBGL"), open="r")
coex <- fromGXL(con)
close(con)

biConnComp(coex)
articulationPoints(coex)
}
\keyword{ models }
