\name{tsort}
\alias{tsort}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ topological sort of vertices of a digraph }
\description{
returns vector of zero-based indices of vertices
of a DAG in topological sort order
}
\usage{
tsort(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ instance of class graf (definition supplied with RBGL)}
}
\details{
calls to the topological\_sort algorithm of BGL.  will check in BGL
whether the input is a DAG and will report to cerr and return
zeroes if it is not.
}
\value{
a numerical vector enumerating vertices in the
topological sort sequence, 0-based
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{Very preliminary}

\examples{
data(FileDep)
vertices(FileDep)[tsort(FileDep)+1]
FD2 <- FileDep
# now introduce a cycle
FD2@edges[[20]] <- c("bar_cpp", "dax_h")
tsort(FD2)
}
\keyword{ models }% at least one, from doc/KEYWORDS
