\name{highlyConnSG}
\alias{highlyConnSG}
\title{Compute highly connected subgraphs for an undirected graph}
\description{Compute highly connected subgraphs for an undirected graph}
}
\usage{
highlyConnSG(g, sat=3, ldv=c(3,2,1))
}
\arguments{
  \item{g}{an instance of the \code{graph} class with \code{edgemode}
    \dQuote{undirected}}
  \item{sat}{singleton adoption threshold, positive integer }
  \item{ldv}{heuristics to remove lower degree vertice, a decreasing sequence of positive integer }
}
\details{
A graph G with n vertices is highly connected if its connectivity k(G) > n/2.  The HCS algorithm partitions a given graph into a set of highly connected subgraphs, by using minimum-cut algorithm recursively.  To improve performance, the approach is refined by adopting singletons, removing low degree vertices and merging clusters.    

On singleton adoption: 
after each round of partition,  some highly connected subgraphs could be
singletons (i.e., a subgraph contains only one node).
To reduce the number of singletons, therefore reduce number of clusters, 
we try to get "normal" subgraphs to "adopt" them.  If a singleton, s, has n 
neighbours in a highly connected subgraph c, and n > sat, we add s to c.  
To adapt to the modified subgraphs, this adoption process is repeated until 
no further such adoption. 

On lower degree vertices: when the graph has low degree vertices, minimum-cut
algorithm will just repeatedly separate these vertices from the rest. 
To reduce such expensive and non-informative computation, we "remove" these 
low degree vertices first before applying minimum-cut algorithm.  
Given ldv = (d1, d2, ..., dp), (d[i] > d[i+1] > 0), we repeat the following
(i from 1 to p): remove all the highly-connected-subgraph found so far; 
remove vertices with degrees < di; find highly-connected-subgraphs; 
perform singleton adoptions. 
}
\value{
  A list of clusters, each is given as vertices in the graph.
}
\references{ A Clustering Algorithm based on Graph Connectivity by E. Hartuv, R. Shamir, 1999.  }
\author{Li Long <li.long@isb-sib.ch>}

\seealso{\code{\link{edgeConnectivity}}, \code{\link{minCut}}  }

\examples{
coex <- fromGXL(file(system.file("XML/hcs.gxl",package="RBGL"), open="r"))
highlyConnSG(coex)
}
\keyword{ models }
