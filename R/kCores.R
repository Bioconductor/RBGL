kCores <- function(g, EdgeType=c("in", "out"))
{
   nv <- length(nodes(g))
   core <- array(0, length(nodes(g)), dimnames=list(nodes(g)))

   # compute the degrees of vertices
   # order the set of vertices V in increasing order of their degrees
   if ( edgemode(g) == "undirected")
       deg <- sort(degree(g))
   else   # directed 
   {
        indeg <- core
        outdeg <- core
	if ( any(EdgeType == "in") )  indeg <- degree(g)$inDegree
	if ( any(EdgeType == "out") )  outdeg <- degree(g)$outDegree
	deg <- indeg + outdeg
        deg <- sort(deg)
   }

   # for each v in V in the order do
   for ( i in 1:nv )
   {
       v = names(deg)[i]

       # core[v] = degree[v]
       core[v] = deg[v]

       # for each u in neighbors(v) do
       if ( edgemode(g) == "undirected") 
           ul <- adj(g, v)[[1]]
       else # directed 
       {
           innb <- list()
           outnb <- list()
	   if ( any(EdgeType == "in") )  outnb <- edges(g, v)[[1]]
	   if ( any(EdgeType == "out") )  innb <- inEdges(v, g)[[1]]
           ul <- unlist(c(innb, outnb))
       }

       for ( u in ul )
       {
	  # if degree[u] > degree[v]
          if ( deg[u] > deg[v] ) 
	  {
	     # degree[u] = degree[u] - 1
	     deg[u] = deg[u] - 1

	     # reorder V accordingly
   	     deg <- sort(deg)
	  }
       }
   }
   core
}
