tsort <- function(x) {
 if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 nv <- length(nodes(x))
 ne <- length(unlist(edges(x)))
 .Call("BGL_tsort_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(x)))
}


mstree.kruskal <- function(x) {
 if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 nv <- length(nodes(x))
 ne <- length(unlist(edges(x)))
 is.directed <- as.integer(edgemode(x) == "directed")
 if (is.directed)
 {
 	ans <- .Call("BGL_KMST_D", 
      		as.integer(nv), as.integer(ne),
      		as.integer(adjListBGL(x)), as.double(unlist(edgeWeights(x))))
 }
 else {
 	ans <- .Call("BGL_KMST_U", 
      		as.integer(nv), as.integer(ne),
      		as.integer(adjListBGL(x)), as.double(unlist(edgeWeights(x))))
 }
 names(ans) <- c("edgeList", "weights")
 ans$nodes <- nodes(x)
 ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans
}

bfs <- function(x,init.ind=1) {
 if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 if (init.ind < 1) stop("use 1-based counting for init.ind")
 nv <- length(nodes(x))
 if (init.ind > nv) stop(paste("only",nv,"nodes but init.ind is ",init.ind,sep=" "))
 ne <- length(unlist(edges(x)))
 ans <- .Call("BGL_bfs_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(x)), as.integer(unlist(edgeWeights(x))),
      as.integer(init.ind-1))
# names(ans) <- c("edgeList", "weights")
# ans$nodes <- nodes(x)
# ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans+1
}

require(graph)
if (!isGeneric("dfs")) setGeneric("dfs", function(object)standardGeneric("dfs"))
setMethod("dfs", "graphNEL", function(object) {
 nv <- length(nodes(object))
 ne <- length(unlist(edges(object)))
 ans <- .Call("BGL_dfs_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(object)), as.double(unlist(edgeWeights(object))))
 names(ans) <- c("discovered", "finish")
 lapply(ans,function(x)x+1)
})

setMethod("dfs", "graph", function(object)
 stop("dfs only defined at present for graphNEL"))


dijkstra.sp <- function(x,init.ind=1) {
 if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 if (init.ind < 1) stop("use 1-based counting for init.ind")
 nv <- length(nodes(x))
 if (init.ind > nv) stop(paste("only",nv,"nodes but init.ind is ",init.ind,sep=" "))
 ne <- length(unlist(edges(x)))
 ans <- .Call("BGL_dijkstra_shortest_paths_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(x)), as.double(unlist(edgeWeights(x))),
      as.integer(init.ind-1))
 names(ans) <- c("distances", "penult")
 ans$penult <- ans$penult+1
 ans[["start"]] <- init.ind
 ans
}
