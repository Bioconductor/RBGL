\name{isomorphism}
\alias{isomorphism}

\title{Compute isomorphism from vertices in one graph to those in another graph}

\description{Compute isomorphism from vertices in one graph to those in another graph}

\usage{
isomorphism(g1, g2)
}

\arguments{
  \item{g1}{one instance of the \code{graph} class }
  \item{g2}{one instance of the \code{graph} class }
}

\details{
As stated in documentation on isomorphism in Boost Graph Library:
An isomorphism is a 1-to-1 mapping of the vertices in one graph to the vertices 
of another graph such that adjacency is preserved. Another words, given graphs 
G1 = (V1,E1) and G2 = (V2,E2) an isomorphism is a function f such that for all 
pairs of vertices a,b in V1, edge (a,b) is in E1 if and only if edge (f(a),f(b))
is in E2.

}

\value{
  Output is true if there exists an isomorphism between g1 and g2, otherwise it's false.
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{ }

\examples{
con1 <- file(system.file("XML/dijkex.gxl",package="RBGL"), open="r")
g1 <- fromGXL(con1)
close(con1)

con2 <- file(system.file("XML/conn2.gxl",package="RBGL"), open="r")
g2 <- fromGXL(con2)
close(con2)

isomorphism(g1, g2)
}
\keyword{ models }
