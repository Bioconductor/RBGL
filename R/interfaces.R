tsort <- function(x) {
# if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 nv <- length(nodes(x))
 ne <- length(unlist(edges(x)))
 .Call("BGL_tsort_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(x)))
}

mstree.kruskal <- function(x) {
# if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 nv <- length(nodes(x))
 ne <- length(unlist(edges(x)))
# is.directed <- as.integer(edgemode(x) == "directed")
# if (is.directed)
# {
 	ans <- .Call("BGL_KMST_D", 
      		as.integer(nv), as.integer(ne),
      		as.integer(adjListBGL(x)), as.double(unlist(edgeWeights(x))))
# }
# else {
#
# here we deal with the undirected graph representation, for
# boost the adjacency list DOES NOT have reciprocal entries
# if it has recip entries you should use a directed representation
# for boost 
#
#        adjlist <- adjListBGL(x)
#        wgts <- unlist(edgeWeights(x))
# 	ans <- .Call("BGL_KMST_U", 
#      		as.integer(nv), as.integer(ne),
#      		as.integer(adjlist), as.double(wgts))
# }
 names(ans) <- c("edgeList", "weights")
 ans$nodes <- nodes(x)
 ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans
}

bfs <- function(x,init.ind=1) {
# if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
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

if (!isGeneric("dfs")) setGeneric("dfs", function(object)standardGeneric("dfs"))
setMethod("dfs", "graph", function(object) {
 nv <- length(nodes(object))
 ne <- length(unlist(edges(object)))
 ans <- .Call("BGL_dfs_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(object)), as.double(unlist(edgeWeights(object))))
 names(ans) <- c("discovered", "finish")
 lapply(ans,function(x)x+1)
})

#setMethod("dfs", "graph", function(object)
# stop("dfs only defined at present for graphNEL"))


dijkstra.sp <- function(x,init.ind=1) {
#    if (class(x) != "graphNEL") 
#        stop("presently only works for class graphNEL from Bioconductor graph library")
    init.ind.ok <- init.ind
    if (is.character(init.ind)) 
        if (init.ind %in% nodes(x)) 
            init.ind.ok <- (1:length(nodes(x)))[init.ind == nodes(x)]
        else stop("character init.ind not in nodes(x)")
    else if (init.ind < 1) 
        stop("use 1-based counting for init.ind")
    nv <- length(nodes(x))
    if (init.ind.ok > nv) 
        stop(paste("only", nv, "nodes but init.ind is ", init.ind.ok, 
            sep = " "))
    ne <- length(unlist(edges(x)))
    ans <- .Call("BGL_dijkstra_shortest_paths_D", as.integer(nv), 
        as.integer(ne), as.integer(adjListBGL(x)), as.double(unlist(edgeWeights(x))), 
        as.integer(init.ind.ok - 1))
    names(ans) <- c("distances", "penult")
    names(ans[[1]]) <- nodes(x)
    ans$penult <- ans$penult + 1
    ans[["start"]] <- init.ind.ok
    ans
}

sp.between <- function (g, start, finish) 
{
# (c) 2003 VJ Carey, all rights reserved
#
# function uses BGL dijkstra shortest paths
# given s=start node, f=end node in graph g,
# return list with length of shortest path joining
# s and f, and vector giving trajectory
#
# debugged 24 sep03, did not need to recompute
# distance, and did not correctly step through
# penultimates!
#
    f <- finish
    s <- start
    no <- nodes(g)
    if (any(is.na(lk <- match(c(s, f), no)))) 
        stop(paste(paste(c(s, f)[is.na(lk)], collapse = " "), 
            "not in nodes of g"))
    s <- (1:length(no))[no == s]
    f <- (1:length(no))[no == f]
    ff <- f
    sp <- dijkstra.sp(g, s)
    if (sp$distances[ff] >= .Machine$double.xmax)
		stop(paste("no path from",no[s],"to",no[f]))
    pens <- sp$penult
    path <- f
    while (path[1] != s) {
        path <- c(pens[f], path)
        f <- pens[f]
    }
    list(length = sp$distances[ff], path = no[path])
}

connComp <- function(g) 
   {
   if (edgemode(g) == "directed")
    .Call("BGL_connected_components_U", as.integer(nv), 
        as.integer(ne), as.integer(adjListBGL(x)), 
        as.double(unlist(edgeWeights(x))))
   else stop("only supporting undirected presently")
   }
