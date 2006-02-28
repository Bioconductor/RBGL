kCores <- function(g)
{
   nv <- length(nodes(g))
   core <- array(0, length(nodes(g)), dimnames=list(nodes(g)))

   # compute the degrees of vertices
   # order the set of vertices V in increasing order of their degrees
   deg <- sort(degree(g))

   # for each v in V in the order do
   for ( i in 1:nv )
   {
       v = names(deg)[i]

       # core[v] = degree[v]
       core[v] = deg[v]

       # for each u in neighbors(v) do
       ul <- adj(g, v)[[1]]
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
