\name{layout}
\alias{circle.layout}
\alias{circleLayout}
\alias{kamada.kawai.spring.layout}
\alias{kamadaKawaiSpringLayout}
\alias{fruchtermanReingoldForceDirectedLayout}
\alias{gursoyAtunLayout}
\alias{randomGraphLayout}

\title{Layout an undirected graph in 2D}

\description{Layout an undirected graph in 2D}

\usage{
circleLayout(g, radius=1)
kamadaKawaiSpringLayout( g, edge_or_side=1, es_length=1 )
fruchtermanReingoldForceDirectedLayout(g, width=1, height=1)
randomGraphLayout(g, minX=0, maxX=1, minY=0, maxY=1)
}

\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{radius}{radius of a regular n-polygon}

  \item{edge_or_side}{boolean indicating the length is for an edge or for a 
    side, default is for an edge }
  \item{es_length}{the length of an edge or a side for layout }

  \item{width}{the width of the dislay area, all x coordinates fall in [-width/2, width/2]}
  \item{height}{the height of the display area, all y coordinates fall in [-height/2, height/2]}

  \item{minX}{minimum x coordinate}
  \item{maxX}{maximum x coordinate}
  \item{minY}{minimum y coordinate}
  \item{maxY}{maximum y coordinate}
  
}

\details{
If you want to simply draw a graph, you should consider using package 
\emph{Rgraphviz}.  The layout options in package \emph{Rgraphviz}: \code{neato},
\code{circo} and \code{fdp}, correspond to \code{kamadaKawaiSpringLayout}, 
\code{circleLayout} and \code{fruchtermanReingoldForceDirectedLayout}, 
respectively. 

Function \code{circleLayout} layouts the graph with the vertices at the points 
of a regular n-polygon.  The distance from the center of the polygon to each 
point is determined by the \code{radius} parameter.

Function \code{kamadaKawaiSpringLayout} provides Kamada-Kawai spring layout for 
connected, undirected graphs.  User provides either the unit length e of an 
edge in the layout or the length of a side s of the display area.

Function \code{randomGraphLayout} places the points of the graph at random locations.

Function \code{fruchtermanReingoldForceDirectedLayout} performs layout of 
unweighted, undirected graphs.  It's a force-directed algorithm.  The BGL 
implementation doesn't handle disconnected graphs very well, since it doesn't
explicitly give each connected component a region proportional to its size. 

%Function \code{gursoyAtunLayout} performs layout by distributing vertices 
%within a topology.  It's based on self-organized maps.

See documentation on this function in Boost Graph Library for more details.
}

\value{
  A (2 x n) matrix, where n is the number of nodes in the graph, each column
gives the (x, y)-coordinates for the corresponding node.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{ \code{\link[Rgraphviz]{layoutGraph}} }

\examples{
con <- file(system.file("XML/conn.gxl",package="RBGL"), open="r")
coex <- fromGXL(con)
close(con)

coex <- ugraph(coex)

circleLayout(coex)

kamadaKawaiSpringLayout(coex)

randomGraphLayout(coex)

fruchtermanReingoldForceDirectedLayout(coex, 10, 10)

}
\keyword{ models }
