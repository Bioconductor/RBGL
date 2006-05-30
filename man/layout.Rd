\name{layout}
\alias{circle.layout}
\alias{kamada.kawai.spring.layout}
\alias{fruchtermanReingoldForceDirectedLayout}
\alias{gursoyAtunLayout}
\alias{randomGraphLayout}
\title{Layout an undirected graph in 2D}
\description{Layout an undirected graph in 2D}
}
\usage{
circle.layout(g, radius=1)
kamada.kawai.spring.layout( g, edge_or_side=1, es_length=1 )
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{radius}{radius of a regular n-polygon}
  \item{edge_or_side}{boolean indicating the length is for an edge or for a 
    side, default is for an edge }
  \item{es_length}{the length of an edge or a side for layout }
}
\details{
According to the documentation on these algorithms in Boost Graph Library:
\code{circle.layout} layouts the graph with the vertices at the points of a 
regular n-polygon.  The distance from the center of the polygon to each point 
is determined by the \code{radius} parameter.

\code{kamada.kawai.spring.layout} provides Kamada-Kawai spring layout for 
connected, undirected graphs.  User provides either the unit length e of an 
edge in the layout or the length of a side s of the display area.

\code{randomGraphLayout} places the points of the graph at random locations.

\code{fruchtermanReingoldForceDirectedLayout} performs layout of unweighted, 
undirected graphs.  It handles disconnected graphs.  It's a force-directed 
algorithm. 

\code{gursoyAtunLayout} performs layout by distributing vertices within a 
topology.  It's based on self-organized maps.

See documentation on this function in Boost Graph Library for more details.
}
\value{
  A list of the (x, y)-coordinates for vertices in the graph.
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex <- ugraph(coex)
circle.layout(coex)
kamada.kawai.spring.layout(coex)
randomGraphLayout(coex)
fruchtermanReingoldForceDirectedLayout(coex, 10, 10)
}
\keyword{ models }
