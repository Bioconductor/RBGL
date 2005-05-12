\name{layout}
\alias{circleLayout}
\alias{kamadaKawaiSpringLayout}
\title{Layout an undirected graph in 2D}
\description{Layout an undirected graph in 2D}
}
\usage{
circleLayout(g, radius=1)
kamadaKawaiSpringLayout( g, edge_or_side=1, es_length=1 )
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
\code{circleLayout} layouts the graph with the vertices at the points of a 
regular n-polygon.  The distance from the center of the polygon to each point 
is determined by the \code{radius} parameter.

\code{kamadaKawaiSpringLayout} provides Kamada-Kawai spring layout for 
connected, undirected graphs.  User provides either the unit length e of an 
edge in the layout or the length of a side s of the display area.
}
\value{
  A list of the (x, y)-coordinates for vertices in the graph.
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex@edgemode <- "undirected"
circleLayout(coex)
kamadaKawaiSpringLayout(coex)
}
\keyword{ models }
