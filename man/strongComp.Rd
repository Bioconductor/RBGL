\name{strongComp}
\alias{strongComp}
\title{Identify Strongly Connected Components}

\description{The strongly connected components in
a directed graph are identified and returned as a list.
}
\usage{
strongComp(g)
}
\arguments{
  \item{g}{ graph with \code{edgemode} \dQuote{directed}.}
}
\details{
Tarjan's algorithm is used to determine all strongly connected
components of a \emph{directed graph}.
}
\value{
A list whose length is the number of strongly connected components in
\code{g}. Each element of the list is a vector of the node labels for
the nodes in that component.
}

\references{Boost Graph Library  }
\author{Vince Carey <stvjc@channing.harvard.edu>}

\seealso{\code{\link[graph]{connComp}},\code{\link{connectedComp}}, \code{\link{same.component}} }

\examples{
km <- fromGXL(file(system.file("XML/kmstEx.gxl",package="RBGL")))
km@nodes <- c(km@nodes,"F","G","H")
km@edgeL$F <- list(edges=numeric(0))
km@edgeL$G <- list(edges=8)
km@edgeL$H <- list(edges=7)
strongComp(km)
connectedComp(ugraph(km))
}
\keyword{ models }
