\name{is.triangulated}
\alias{is.triangulated}
\title{Decide if a graph is triangulated}
\description{Decide if a graph is triangulated}
}
\usage{
is.triangulated(g)
}
\arguments{
  \item{g}{an instance of the \code{graph} class }
}
\details{
Our implementation is based on the following:

Combinatorial Optimization: algorithms and complexity  (p. 403)
by C. H. Papadimitriou, K. Steiglitz

A graph G = (V, E) is triangulated (i.e. chordal) if all cycles
[v1, v2, ..., vk] of length 4 or more have a chord, i.e., an edge
[vi, vj] with j != i +/- 1 (mod k)

An equivalent definition of chordal graphs is:

G is chordal iff either G is an empty graph, or
there is an v in V such that
(i) the neighborhood of v (i.e., v and its adjacent nodes) forms a clique, and
(ii) recursively, G-v is chordal

}
\value{
Output is true if \code{g} is triangulated, false otherwise.
}
\references{}
\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
is.triangulated(coex)
coex <- fromGXL(file(system.file("XML/hcs.gxl",package="RBGL")))
is.triangulated(coex)
}
\keyword{ models }
