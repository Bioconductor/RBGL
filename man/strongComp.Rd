\name{strongComp}
\alias{strongComp}

\title{Identify Strongly Connected Components}

\description{The strongly connected components in
a directed graph are identified and returned as a list.  }

\usage{
strongComp(g)
}

\arguments{
  \item{g}{ graph with \code{edgemode} \dQuote{directed}.}
}

\details{
Tarjan's algorithm is used to determine all strongly connected
components of a \emph{directed graph}.
}

\value{
A list whose length is the number of strongly connected components in
\code{g}. Each element of the list is a vector of the node labels for
the nodes in that component.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Vince Carey <stvjc@channing.harvard.edu>}

\seealso{\code{\link[graph]{connComp}},\code{\link{connectedComp}}, \code{\link{same.component}} }

\examples{
con <- file(system.file("XML/kmstEx.gxl",package="RBGL"), open="r")
km <- fromGXL(con)
close(con)

km<- graph::addNode(c("F","G","H"), km)
km<- addEdge("G", "H", km, 1)
km<- addEdge("H", "G", km, 1)
strongComp(km)
connectedComp(ugraph(km))
}
\keyword{ models }
