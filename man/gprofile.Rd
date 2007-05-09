\name{gprofile}
\alias{gprofile}
\title{Compute profile for a graph}

\description{Compute profile for a graph}

\usage{
gprofile(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class} 
}

\details{
The profile of a given graph is the sum of bandwidths for all the vertices in 
the graph.

See documentation on this function in Boost Graph Library for more details.
}

\value{
\item{profile}{the profile of the graph}
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
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL"), open="r"))
gprofile(coex)
}
\keyword{ models }
