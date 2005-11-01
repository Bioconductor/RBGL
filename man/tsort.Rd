\name{tsort}
\alias{tsort}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ topological sort of vertices of a digraph }
\description{
returns vector of zero-based indices of vertices
of a DAG in topological sort order
}
\usage{
tsort(x) # now x assumed to be Bioconductor graph graphNEL
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ instance of class graphNEL from Bioconductor graph class}
}
\details{
calls to the topological\_sort algorithm of BGL.  will check in BGL
whether the input is a DAG and return a vector of
zeroes (of length length(nodes(x))) if it is not.  Thus this function can be used to check for
cycles in a digraph.
}
\value{
a numerical vector enumerating vertices in the
topological sort sequence, 0-based
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\examples{
 data(FileDep)
 tsind <- tsort(FileDep)
 tsind
 nodes(FileDep)[tsind+1]
 FD2 <- FileDep
 # now introduce a cycle
 FD2 <- addEdge("bar_o", "dax_h", FD2, 1)
 tsort(FD2)
 }
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
