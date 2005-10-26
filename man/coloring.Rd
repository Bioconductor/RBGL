\name{sequential.vertex.coloring}
\alias{sequential.vertex.coloring}
\title{Compute how many colors are needed to color a graph}
\description{Compute how many colors are needed to color a graph}
}
\usage{
sequential.vertex.coloring(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class} 
}
\details{
TO BE FILLED IN
}
\value{
\item{profile}{the profile of the graph}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
sequential.vertex.coloring(coex)
}
\keyword{ models }
