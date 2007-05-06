\name{biConnComp}
\alias{biConnComp}
\alias{articulationPoints}
\title{Compute biconnected components for a graph}
\description{Compute biconnected components for a graph}
}
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
\item{no. of biconnected components}{number of biconnected components in the graph}
\item{biconnected components}{component label for each edge in the graph}

For \code{articulationPoints}:
\item{no. of articulation points}{number of articulation points in the graph}
\item{articulation points}{list of articulation points in the graph}

}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL"), open="r"))
biConnComp(coex)
articulationPoints(coex)
}
\keyword{ models }
