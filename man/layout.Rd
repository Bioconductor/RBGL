\name{layout}
\alias{layout}
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
\code{circle.layout} layouts the graph with the vertices at the points of a regular n-polygon.  The distance from the center of the polygon to each point is determined by the \code{radius} parameter.

\code{kamada.kawai.spring.layout} provides Kamada-Kawai spring layout for connected, undirected graphs.  User provides either the unit length e of an edge in the layout or the length of a side s of the display area.
}
\value{
\item{circle.layout}{ the (x, y)-coordinates of the vertices }
\item{kamada.kawai.spring.layout}{ the (x, y)-coordinates of the vertices }
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
coex@edgemode <- "undirected"
circle.layout(coex)
kamada.kawai.spring.layout(coex)
}
\keyword{ models }
