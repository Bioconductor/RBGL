\name{johnson.all.pairs.sp}
\alias{johnson.all.pairs.sp}
\title{compute shortest path distance matrix for all pairs of nodes}
\description{
compute shortest path distance matrix for all pairs of nodes
}
\usage{
johnson.all.pairs.sp(g)
}
\arguments{
  \item{g}{ graph object for which edgeMatrix and edgeWeights are defined}
}
\details{
Uses BGL algorithm.
}
\value{
matrix of shortest path lengths, read from row node to col node
}
\references{BGL  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{  }


\seealso{  }

\examples{
z <- fromGXL(file(system.file("dot/joh.gxl", package="RBGL")))
johnson.all.pairs.sp(z)
}
\keyword{ models }
