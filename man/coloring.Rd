\name{sequential.vertex.coloring}
\alias{sequential.vertex.coloring}
\title{Compute a vertex coloring for a graph}
\description{Compute vertex coloring for a graph}
}
\usage{
sequential.vertex.coloring(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class} 
}
\details{
A vertex coloring for a graph is to assign a color for each vertex so that no 
two adjacent vertices are of the same color.  We designate the colors as 
sequential integers: 1, 2, ....

For ordered vertices, \code{v1}, \code{v2}, ..., \code{vn}, for k = 1, 2, ..., 
n, this algorithm assigns \code{vk} to the smallest possible color.  It does 
NOT guarantee to use minimum number of colors.

See documentations on these algorithms in Boost Graph Library for more details.
}
\value{
\item{no. of colors needed}{how many colors to use to color the graph}
\item{colors of nodes} { color label for each vertex }
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
sequential.vertex.coloring(coex)
}
\keyword{ models }
