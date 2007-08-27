\name{removeSelfLoops}
\alias{removeSelfLoops}

\title{remove self loops in a graph}

\description{remove self loops in a graph}

\usage{
removeSelfLoops(g)
}

\arguments{
  \item{g}{one instance of the \code{graph} class }
}

\details{
If a given graph contains self-loop(s), \code{removeSelfLoops} removes them.
This is for those functions that cannot handle graphs with self-loops.
}

\value{
  A new graph without self loops.
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{ }

\examples{
con <- file(system.file("XML/dijkex.gxl",package="RBGL"))
g1 <- fromGXL(con)
close(con)

g2 <- ugraph(g1)
removeSelfLoops(g2)
}
\keyword{ models }
