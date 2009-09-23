\name{RBGL.overview}
\alias{RBGL.overview}
\title{RBGL.overview}

\description{The RBGL package consists of a number of interfaces to the Boost
C++ library for graph algorithms.  This page follows, approximately,
the chapter structure of the monograph on the Boost Graph Library
by Siek et al., and gives hyperlinks to documentation on
R functions currently available, along with the names
of formal parameters to these functions. }

\section{ basicAlgs }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL]{bandwidth}} \tab \code{g}\cr
\code{\link[RBGL]{bfs}} \tab \code{object,node,checkConn}\cr
\code{\link[RBGL:bfs]{dfs}} \tab \code{object,node,checkConn}\cr
\code{\link[RBGL:edgeConn]{edgeConnectivity}} \tab \code{g}\cr
\code{\link[RBGL]{gprofile}} \tab \code{g}\cr
\code{\link[RBGL]{isomorphism}} \tab \code{g1,g2}\cr
\code{\link[RBGL]{minCut}} \tab \code{g}\cr
\code{\link[RBGL:transClosure]{transitive.closure}} \tab \code{g}\cr
\code{\link[RBGL]{tsort}} \tab \code{x}\cr
}
}

\section{ ShortestPaths }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL]{bellman.ford.sp}} \tab \code{g,start}\cr
\code{\link[RBGL]{dag.sp}} \tab \code{g,start}\cr
\code{\link[RBGL]{dijkstra.sp}} \tab \code{g,start}\cr
\code{\link[RBGL]{extractPath}} \tab \code{s,f,pens}\cr
\code{\link[RBGL]{johnson.all.pairs.sp}} \tab \code{g}\cr
\code{\link[RBGL]{sp.between}} \tab \code{g,start,finish}\cr
\code{\link[RBGL:sp.between]{sp.between.old}} \tab \code{g,start,finish}\cr
\code{\link[RBGL:sp.between]{sp.between.scalar}} \tab \code{g,start,finish}\cr
}
}
\section{ MinimumSpanningTree }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL]{mstree.kruskal}} \tab \code{x}\cr
\code{\link[RBGL:mstree.prim]{prim.minST}} \tab \code{g}\cr
}
}
\section{ ConnectedComponents }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL]{connectedComp}} \tab \code{g}\cr
\code{\link[RBGL]{highlyConnSG}} \tab \code{g,sat,ldv}\cr
\code{\link[RBGL:incrConnComp]{incremental.components}} \tab \code{g}\cr
\code{\link[RBGL:incrConnComp]{init.incremental.components}} \tab \code{g}\cr
\code{\link[RBGL:incrConnComp]{same.component}} \tab \code{g,node1,node2}\cr
\code{\link[RBGL]{strongComp}} \tab \code{g}\cr
}
}
\section{ MaximumFlow }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL:maxFlow]{edmonds.karp.max.flow}} \tab \code{g,source,sink}\cr
\code{\link[RBGL:maxFlow]{push.relabel.max.flow}} \tab \code{g,source,sink}\cr
}
}
\section{ SparseMatrixOrdering }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL:ordering]{cuthill.mckee.ordering}} \tab \code{g}\cr
\code{\link[RBGL:ordering]{minDegreeOrdering}} \tab \code{g,delta}\cr
\code{\link[RBGL:ordering]{sloan.ordering}} \tab \code{g,w1,w2}\cr
}
}
\section{ LayoutAlgorithms }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL:layout]{circle.layout}} \tab \code{g,radius}\cr
\code{\link[RBGL:layout]{kamada.kawai.spring.layout}} \tab \code{g,edge_or_side,es_length}\cr
}
}
\section{ GraphClustering }{
\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL:bccluster]{betweenness.centrality.clustering}} \tab \code{g,threshold,normalize}\cr
}
}

\section{ Betweenness }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL:betweenness]{brandes.betweenness.centrality}} \tab \code{g}\cr
}
}
\section{ Wavefront }{

\tabular{ll}{
Functions \tab parameters \cr 
\code{\link[RBGL:wavefront]{aver.wavefront}} \tab \code{g}\cr
\code{\link[RBGL:wavefront]{ith.wavefront}} \tab \code{g,start}\cr
\code{\link[RBGL:wavefront]{maxWavefront}} \tab \code{g}\cr
\code{\link[RBGL:wavefront]{rms.wavefront}} \tab \code{g}\cr
}
}

\references{
Boost Graph Library ( www.boost.org/libs/graph/doc/index.html )

The Boost Graph Library: User Guide and Reference Manual;
by Jeremy G. Siek, Lie-Quan Lee, and Andrew Lumsdaine;
(Addison-Wesley, Pearson Education Inc., 2002), xxiv+321pp.
ISBN 0-201-72914-8
}

\keyword{models}

