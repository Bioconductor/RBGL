\name{astarSearch}
\alias{astarSearch}
\title{Compute astarSearch for a graph}
\description{Compute astarSearch for a graph}
}
\usage{
astarSearch(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class} 
}
\details{
NOT IMPLEMENTED YET.  TO BE FILLED IN
}
\value{
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL"), open="r"))
astarSearch(coex)
}
\keyword{ models }
