\name{extractPath}
\alias{extractPath}
\title{convert a dijkstra.sp predecessor
structure into the path joining two nodes   }
\description{
determine a path between two nodes in a graph,
using output of \code{\link{dijkstra.sp}}.
}
\usage{
extractPath(s, f, pens)
}
\arguments{
  \item{s}{ s: index of starting node in nodes vector
of the graph from which \code{pens} was derived}
  \item{f}{ f: index of ending node in nodes vector }
  \item{pens}{ pens: predecessor index vector as returned
in the \code{preds} component of \code{\link{dijkstra.sp}} output}
}
%\details{
%}
%\value{
%}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


\seealso{  \code{\link[e1071]{allShortestPaths}}}

\examples{
data(FileDep)
dd <- dijkstra.sp(FileDep)
extractPath(1,9,dd$pen)
}
\keyword{ models }
