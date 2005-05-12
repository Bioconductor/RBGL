\name{extractPath}
\alias{extractPath}
\title{convert a dijkstraSP predecessor structure into the path joining two nodes   }
\description{
determine a path between two nodes in a graph,
using output of \code{\link{dijkstraSP}}.
}
\usage{
extractPath(s, f, pens)
}
\arguments{
  \item{s}{ index of starting node in nodes vector
of the graph from which \code{pens} was derived}
  \item{f}{ index of ending node in nodes vector }
  \item{pens}{ predecessor index vector as returned
in the \code{preds} component of \code{\link{dijkstraSP}} output}
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
dd <- dijkstraSP(FileDep)
extractPath(1,9,dd$pen)
}
\keyword{ models }
