\name{FileDep}
\alias{FileDep}
\non_function{}
\title{ FileDep: a graphNEL object representing a file dependency dataset
example in boost graph library }

\description{ FileDep: a graphNEL object representing a file dependency dataset
example in boost graph library }

\usage{data(FileDep)}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\value{an instance of graphNEL}

\examples{
# this is how the graph of data(FileDep) was obtained
library(graph)
fd <- file(system.file("XML/FileDep.gxl",package="RBGL"), open="r")
show(fromGXL(fd))
if (require(Rgraphviz))
 {
 data(FileDep)
 plot(FileDep)
 }
close(fd)
}
\keyword{graphs}
