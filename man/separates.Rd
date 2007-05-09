\name{separates}
\alias{separates}

\title{ A function to test whether a subset of nodes separates
 two other subsets of nodes. }

\description{
The function tests to see whether a set of nodes, \code{S1}, separates 
all nodes in \code{a} from all nodes in \code{b}.
}

\usage{
separates(a, b, S1, g)
}

\arguments{
  \item{a}{ The names of the nodes in the from set. }
  \item{b}{ The names of the nodes in the to set. }
  \item{S1}{ The names of the nodes in the separation set. }
  \item{g}{ An instance of the \code{graph} class. All nodes
     named in the other arguments must be nodes of this graph. }
}

\details{
 The algorithm is quite simple. A subgraph is created by removing the
 nodes named in \code{S1} from \code{g}. Then all paths between
 elements of \code{a} to elements of \code{b} are tested for. If any
 path exists the function returns \code{FALSE}, otherwise it returns 
\code{TRUE}.
}

\value{
 Either \code{TRUE} or \code{FALSE} depending on whether \code{S1} separates
\code{a} from \code{b} in \code{g1}.
}

\references{ S. Lauritzen, Graphical Models, OUP. }

\author{ R. Gentleman }

\seealso{ \code{\link{johnson.all.pairs.sp}}}
\examples{
  km <- fromGXL(file(system.file("XML/kmstEx.gxl",package="RBGL")))
  separates("B", "A", "E", km)
  separates("B", "A", "C", km)
}
\keyword{manip }
