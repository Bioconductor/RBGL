\name{tsort}
\alias{tsort}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{ topological sort of vertices of a digraph }

\description{ returns vector of zero-based indices of vertices
of a DAG in topological sort order }

\usage{
tsort(x) # now x assumed to be Bioconductor graph graphNEL
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{ instance of class graphNEL from Bioconductor graph class}
}

\details{
calls to the topological\_sort algorithm of BGL.  will check in BGL
whether the input is a DAG and return a vector of zeroes (of length 
length(nodes(x))) if it is not.  Thus this function can be used to check for
cycles in a digraph.
}

\value{
A character vector of vertices in the topological sort sequence.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ VJ Carey <stvjc@channing.harvard.edu> }

\examples{
 data(FileDep)
 tsind <- tsort(FileDep)
 tsind
 FD2 <- FileDep
 # now introduce a cycle
 FD2 <- addEdge("bar_o", "dax_h", FD2, 1)
 tsort(FD2)
 }
\keyword{ graphs }% at least one, from doc/KEYWORDS
