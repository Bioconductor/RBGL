\name{kCliques}
\alias{kCliques}
\title{Find all the k-cliques in an undirected graph}

\description{Find all the k-cliques in an undirected graph }

\usage{
kCliques(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
}

\details{
Notice that there are different definitions of k-clique in different context.

In computer science, a k-clique of a graph is a clique, i.e., a complete subgraph, of k nodes. 

In Social Network Analysis, a k-clique in a graph is a subgraph where the distance between any two nodes is no greater than k.

Here we take the definition in Social Network Analysis.

Let D be a matrix, D[i][j] is the shortest path from node i to node j.  Algorithm is outlined as following: 
(1) use Johnson's algorithm to fill D; let N = max(D[i][j]) for all i, j; 
(2) each edge is a 1-clique by itself; 
(3) for k = 2, ..., N, try to expand each (k-1)-clique to k-clique: 
    (3.1) consider a (k-1)-clique the current k-clique KC;
    (3.2) repeat the following:
          if for all nodes j in KC, D[v][j] <= k, add node v to KC;
    (3.3) eliminate duplicates;
(4) the whole graph is N-clique.
}

\value{
A list of length N; k-th entry (k = 1, ..., N) is a list of all the k-cliques in graph \code{g}. 
}

\references{
Social Network Analysis: Methods and Applications.  By S. Wasserman and K. Faust, pp. 258. 
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
coex <- fromGXL(file(system.file("XML/snacliqueex.gxl",package="RBGL")))
kCliques(coex)
}

\keyword{ models }
