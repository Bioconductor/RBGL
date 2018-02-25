\name{boyerMyrvoldPlanarityTest}
\alias{boyerMyrvoldPlanarityTest}

\title{boyerMyrvoldPlanarityTest}

\description{boyerMyrvoldPlanarityTest description }

\usage{
boyerMyrvoldPlanarityTest(g) 
}

\arguments{
  \item{g}{ instance of class graphNEL from Bioconductor graph class}
}

%\details{
%}

\value{
logical, TRUE if test succeeds
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\author{ Li Long <li.long@isb-sib.ch> }

\examples{
con <- file(system.file("XML/dijkex.gxl",package="RBGL"), open="r")
coex <- fromGXL(con)
edgemode(coex) = "undirected"
boyerMyrvoldPlanarityTest(coex)  # only shows runnability, need better case
}

\keyword{ graphs }

