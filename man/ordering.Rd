\name{Ordering}
\alias{Ordering}
\title{Compute vertex ordering for an undirected graph}
\description{Compute vertex ordering for an undirected graph}
}
\usage{
cuthill.mckeee.ordering(g)
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
  \item{W1} {Heuristical weight for the Sloan algorithm }
  \item{W2} {Heuristical weight for the Sloan algorithm }
}
\details{
The goal of the Cuthill-Mckee (and reverse Cuthill-Mckee) ordering algorithm is
to reduce the bandwidth of a graph by reordering the indices assigned to each
vertex. The Cuthill-Mckee ordering algorithm works by a local minimization of
the i-th bandwidths. The vertices are basically assigned a breadth-first search
order, except that at each step, the adjacent vertices are placed in the queue
in order of increasing degree.
The minimum degree ordering algorithm is a fill-in reduction matrix reordering
algorithm. When solving a system of equations such as A x = b using a sparse
version of Cholesky factorization (which is a variant of Gaussian elimination
for symmetric matrices), often times elements in the matrix that were once zero
become non-zero due to the elimination process. This is what is referred to as
"fill-in", and fill-in is bad because it makes the matrix less sparse and
therefore consume more time and space in later stages of the elimintation and
in computations that use the resulting factorization. Now it turns out that
reordering the rows of a matrix can have a dramatic affect on the amount of
fill-in that occurs. So instead of solving A x = b, we instead solve the
reordered (but equivalent) system (P A PT)(P x) = P b. Finding the optimal
ordering is an NP-complete problem (e.i., it can not be solved in a reasonable
amount of time) so instead we just find an ordering that is "good enough" using
heuristics. The minimum degree ordering algorithm is one such approach.
A symmetric matrix A is typically represented with an undirected graph, however
for this function we require the input to be a directed graph. Each nonzero
matrix entry A(i, j) must be represented by two directed edges (e(i,j) and
e(j,i)) in G.
The goal of the Sloan ordering algorithm is to reduce the profile and the
wavefront of a graph by reordering the indices assigned to each vertex. The
Sloan algorithm needs a start and an end vertex. These vertices can be asigned
manually. But there is also  an algorithm sloan-starting-nodes that provides
usually quite good start and end vertices. Each vertex is asigned with a
priority. This priority is a weighted sum of the distance of the vector to the
end vertex (a global criterium) and is called the current degree of vertex.
This current degree basically reflects the status of the renumbering in the
neighborhood of a vertex (a local criterium). Therefore the Sloan algorithm (in
contrast to-McKee) takes into account local as well as global criteria for the
renumbering sequence. One can play around with the relative weights, but the
default values proposed by Sloan (weight1/weight2=1/2) turn out to be pretty
good in most cases.

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
