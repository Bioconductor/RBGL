\name{wavefront}
\alias{ithWavefront}
\alias{maxWavefront}
\alias{averWavefront}
\alias{rmsWavefront}
\title{Compute the i-th/max/average/rms wavefront for a graph}
\description{Compute the i-th/max/average/rms wavefront for a graph}
}
\usage{
ithWavefront(g, start)
maxWavefront(g)
averWavefront(g)
rmsWavefront(g)
}
\arguments{
  \item{start}{a vertex of the \code{graph} class }
  \item{g}{an instance of the \code{graph} class }
}
\details{
Assorted functions on wavefront of a graph.
}
\value{
\item{ithWavefront}{wavefront of the given vertex}
\item{maxWavefront}{maximum wavefront of a graph}
\item{averWavefront}{average wavefront of a graph}
\item{rmsWavefront}{root mean square of all wavefronts}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{edgeConnectivity}}  }

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
ss <- 1
ithWavefront(coex, ss)
maxWavefront(coex)
averWavefront(coex)
rmsWavefront(coex)
}
\keyword{ models }
