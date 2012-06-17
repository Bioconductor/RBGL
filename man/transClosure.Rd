\name{transitive.closure}
\alias{transitive.closure}

\title{Compute transitive closure of a directed graph -- suspended 16 june 2012}

\description{Compute  transitive closure of a directed graph}

\usage{
transitive.closure(g) # not compiling with boost 1.49 -- needs work
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
}

\details{
This function calculates the transitive closure of a directed graph.  
See documentation on this function in Boost Graph Library for more details.
}

\value{
  An object of class \code{graphNEL}.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
con <- file(system.file("XML/dijkex.gxl",package="RBGL"))
coex <- fromGXL(con)
close(con)

transitive.closure(coex)
}
\keyword{ models }
