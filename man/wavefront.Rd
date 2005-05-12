\name{wavefront}
\alias{ith.wavefront}
\alias{maxWavefront}
\alias{aver.wavefront}
\alias{rms.wavefront}
\title{Compute the i-th/max/average/rms wavefront for a graph}
\description{Compute the i-th/max/average/rms wavefront for a graph}
}
\usage{
ith.wavefront(g, start)
maxWavefront(g)
aver.wavefront(g)
rms.wavefront(g)
}
\arguments{
  \item{start}{a vertex of the \code{graph} class }
  \item{g}{an instance of the \code{graph} class }
}
\details{
Assorted functions on wavefront of a graph.
}
\value{
\item{ith.wavefront}{wavefront of the given vertex}
\item{maxWavefront}{maximum wavefront of a graph}
\item{aver.wavefront}{average wavefront of a graph}
\item{rms.wavefront}{root mean square of all wavefronts}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{edgeConnectivity}}  }

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
ss <- 1
ith.wavefront(coex, ss)
maxWavefront(coex)
aver.wavefront(coex)
rms.wavefront(coex)
}
\keyword{ models }
