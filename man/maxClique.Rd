\name{maxClique}
\alias{maxClique}

\title{Find all the cliques in a graph}

\description{Find all the cliques in a graph }

\usage{
maxClique(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
}

\details{
Notice the maximum clique problem is NP-complete, which means it 
cannot be solved by any known polynomial algorithm.

We implemented the algorithm by C. Bron and J. Kerbosch, 
}

\value{
\item{maxClique}{list of all cliques in \code{g}}
}

\references{
Finding all cliques of an undirected graph, 
by C. Bron and J. Kerbosch, 
Communication of ACM, Sept 1973, Vol 16, No. 9.
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
con1 <- file(system.file("XML/conn.gxl",package="RBGL"), open="r")
coex <- fromGXL(con1)
close(con1)

maxClique(coex)

con2 <- file(system.file("XML/hcs.gxl",package="RBGL"), open="r")
coex <- fromGXL(con2)
close(con2)

maxClique(coex)
}
\keyword{ models }
