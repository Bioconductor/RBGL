\name{lambdaSets}
\alias{lambdaSets}

\title{Find all the lambda-sets in an undirected graph}

\description{Find all the lambda-sets in an undirected graph }

\usage{
lambdaSets(g)
}

\arguments{
  \item{g}{an instance of the \code{graph} class }
}

\details{
From reference (1), p. 270:
A set of nodes is a lambda-set if any pair of nodes in the lambda set has 
larger edge connectivity than any pair of nodes consisting of one node from 
within the lamda set and a second node from outside the lamda set.

As stated in reference (2), a lambda set is a maximal subset of nodes who have 
more edge-independent paths connecting them to each other than to outsiders.

A lambda set could be characterized by the minimum edge connectivity \code{k}
among its members, and could be called \code{lambda-k} sets.

Let N be maximum edge connectivity of graph \code{g},
we output all the lambda-k set for all k = 1, ..., N.
}

\value{
Maximum edge connectivity, \code{N}, of the graph \code{g}, and
A list of length N; k-th entry (k = 1, ..., N) is a list of all the lambda-k 
sets in graph \code{g}. 
}

\references{
(1) Social Network Analysis: Methods and Applications.  By S. Wasserman and K. Faust, pp. 269. 
(2) LS sets, lambda sets and other cohesive subsets.  By S. P. Borgatti, M. G. Everett, P. R. Shirey, Social Networks 12 (1990) p. 337-357
}

\author{Li Long <li.long@isb-sib.ch>}

\seealso{}

\examples{
con <- file(system.file("XML/snalambdaex.gxl",package="RBGL"))
coex <- fromGXL(con)
close(con)

lambdaSets(coex)
}

\keyword{ models }
