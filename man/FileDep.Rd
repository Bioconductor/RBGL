\name{FileDep}
\alias{FileDep}
\non_function{}
\title{ FileDep: a graphNEL object representing a file dependency dataset
example in boost graph library }
\description{ FileDep: a graphNEL object representing a file dependency dataset
example in boost graph library }
\usage{#data(FileDep)}
\source{
Boost Graph Library 
}
\references{
}
\examples{
# this is how the graph of data(FileDep) was obtained
library(graph)
fd <- file(system.file("XML/FileDep.gxl",package="RBGL"))
show(fromGXL(fd))
if (require(Rgraphviz))
 {
 data(FileDep)
 plot(FileDep)
 }
}
\keyword{graphs}
