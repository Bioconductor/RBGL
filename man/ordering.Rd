\name{Ordering}
\alias{Ordering}
\title{Compute vertex ordering for an undirected graph}
\description{Compute vertex ordering for an undirected graph}
}
\usage{
cuthill.mckee.ordering(g)
min.degree.ordering(g, delta=0)
sloan.ordering(g, w1=1, w2=2)
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{delta}{Multiple elimination control variable. If it is larger than or 
    equal to zero then multiple elimination is enabled. The value of delta 
    specifies the difference between the minimum degree and the degree of 
    vertices that are to be eliminated.}
  \item{w1} {Heuristical weight for the Sloan algorithm }
  \item{w2} {Heuristical weight for the Sloan algorithm }
}
\details{
According to the documentations on these algorithms in Boost Graph Library:
The goal of the Cuthill-Mckee (and reverse Cuthill-Mckee) ordering algorithm is
to reduce the bandwidth of a graph by reordering the indices assigned to each
vertex. 
The minimum degree ordering algorithm is a fill-in reduction matrix reordering
algorithm. 
The goal of the Sloan ordering algorithm is to reduce the profile and the
wavefront of a graph by reordering the indices assigned to each vertex. 
See documentations on these algorithms in Boost Graph Library for more details.
}
\value{
\item{cuthill.mckee.ordering}{the vertices in the new ordering}
\item{min.degree.ordering}{the vertices in the new ordering}
\item{sloan.ordering}{the vertices in the new ordering}
}
\references{Boost Graph Library by Siek et al.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{ }

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL")))
coex@edgemode <- "undirected"
cuthill.mckee.ordering(coex)
min.degree.ordering(coex)
sloan.ordering(coex)
}
\keyword{ models }
