kCores <- function(g, EdgeType=c("in", "out"))
{
   nv <- numNodes(g)
   core <- array(0, nv, dimnames=list(nodes(g)))

   # compute the degrees of vertices
   # order the set of vertices V in increasing order of their degrees
   if ( edgemode(g) == "undirected")
       deg <- sort(degree(g))
   else   
   {
        indeg <- core
        outdeg <- core
	if ( any(EdgeType == "in") )  indeg <- degree(g)$inDegree
	if ( any(EdgeType == "out") )  outdeg <- degree(g)$outDegree
	deg <- indeg + outdeg
        deg <- sort(deg)
   }

   # for each v in V in the order do
   for ( i in seq_len(nv) )
   {
       v = names(deg)[i]

       core[v] = deg[v]

       # for each u in neighbors(v) do
       if ( edgemode(g) == "undirected") 
           ul <- adj(g, v)[[1]]
       else 
       {
           innb <- list()
           outnb <- list()
	   if ( any(EdgeType == "in") )  outnb <- edges(g, v)[[1]]
	   if ( any(EdgeType == "out") )  innb <- inEdges(v, g)[[1]]
           ul <- unlist(c(innb, outnb))
       }

       for ( u in ul )
       {
          if ( deg[u] > deg[v] ) 
	  {
	     deg[u] = deg[u] - 1
   	     deg <- sort(deg)
	  }
       }
   }
   core
}
