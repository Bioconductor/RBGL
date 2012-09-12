\name{maxClique}
\alias{maxClique}

\title{Find all the cliques in a graph}

\description{Find all the cliques in a graph }

\usage{
maxClique(g, nodes=NULL, edgeMat=NULL)
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
 \item{nodes}{vector of node names, to be supplied if \code{g} is not}
 \item{edgeMat}{2 x p matrix with indices of edges in nodes, one-based,
   only to be supplied if code{g} is not}
}

\details{
Notice the maximum clique problem is NP-complete, which means it 
cannot be solved by any known polynomial algorithm.

We implemented the algorithm by C. Bron and J. Kerbosch, 

It is an error to supply both \code{g} and either of the other arguments.

If \code{g} is not supplied, no checking of the consistency of
\code{nodes} and \code{edgeMat} is performed.
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

%\seealso{}

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
