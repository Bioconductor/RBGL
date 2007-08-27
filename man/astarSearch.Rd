\name{astarSearch}
\alias{astarSearch}
\title{Compute astarSearch for a graph}
\description{Compute astarSearch for a graph}
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
\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual; 
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine; 
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp. 
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
con <- file(system.file("XML/dijkex.gxl",package="RBGL"), open="r")
coex <- fromGXL(con)
close(con)
astarSearch(coex)
}
\keyword{ models }
