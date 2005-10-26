\name{biConnComp}
\alias{biConnComp}
\alias{articulationPoints}
\title{Compute biconnected components for a graph}
\description{Compute biconnected components for a graph}
}
\usage{
biConnComp(g)
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
biConnComp(coex)
}
\keyword{ models }
