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
An undirected graph G = (V, E) is triangulated (i.e. chordal) if all cycles
[v1, v2, ..., vk] of length 4 or more have a chord, i.e., an edge
[vi, vj] with j != i +/- 1 (mod k)

An equivalent definition of chordal graphs is:

G is chordal iff either G is an empty graph, or
there is an v in V such that
\enumerate{
\item the neighborhood of v (i.e., v and its adjacent nodes) forms a clique, and
\item recursively, G-v is chordal
}

}
\value{
The return value is \code{TRUE} if \code{g} is triangulated
 and \code{FALSE} otherwise. An error is thrown if the graph is not 
undirected; you might use \code{\link{ugraph}} to compute the underlying 
graph.
}
\references{
Combinatorial Optimization: algorithms and complexity  (p. 403)
by C. H. Papadimitriou, K. Steiglitz}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/conn.gxl",package="RBGL")))
is.triangulated(coex)
coex <- fromGXL(file(system.file("XML/hcs.gxl",package="RBGL")))
is.triangulated(coex)
}
\keyword{ models }
