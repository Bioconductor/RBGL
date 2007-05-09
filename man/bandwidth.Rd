\name{bandwidth}
\alias{bandwidth}
\title{Compute bandwidth for an undirected graph}
\description{Compute bandwidth for an undirected graph}

\usage{
bandwidth(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
}

\details{
The bandwidth of an undirected graph G=(V, E) is the maximum distance between 
two adjacent vertices. 
See documentation on bandwidth in Boost Graph Library for more details.
}

\value{
\item{bandwidth}{the bandwidth of the given graph}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{  }

\examples{
coex <- fromGXL(file(system.file("XML/dijkex.gxl",package="RBGL"), open="r"))
coex <- ugraph(coex)
bandwidth(coex)
}
\keyword{ models }
