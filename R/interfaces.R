tsort <- function(x) {
 if (edgemode(x) != "directed") stop("requires directed graph")
 nv <- length(nodes(x))
 em <- edgeMatrix(x)
 ne <- ncol(em)
 .Call("BGL_tsort_D", as.integer(nv), as.integer(ne),
      as.integer(em-1))
}

mstree.kruskal <- function(x) {
 nv <- length(nodes(x))
 em <- edgeMatrix(x, duplicates=TRUE) # conform with edgeWeights unlisted
 ne <- ncol(em)
ans <- .Call("BGL_KMST_D", 
      		as.integer(nv), as.integer(ne),
      		as.integer(em-1), as.double(unlist(edgeWeights(x))))
 names(ans) <- c("edgeList", "weights")
 ans$nodes <- nodes(x)
 ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans
}

if (!isGeneric("bfs")) setGeneric("bfs",
        function( graph, init.node=1, checkConn=FALSE) standardGeneric("bfs"))

setMethod("bfs",c("graph", "ANY", "ANY"), 
  function( graph, init.node, checkConn) {
 if (init.node < 1) stop("use 1-based counting for init.node")
 if (checkConn)
   {
   if (length(connectedComp(graph))>1) stop("graph is not connected")
   }
 nv <- length(nodes(graph))
 if (init.node > nv) stop(paste("only",nv,"nodes but init.node is ",init.node,sep=" "))
 em <- edgeMatrix(graph)
 ne <- ncol(em)
 ans <- .Call("BGL_bfs_D", as.integer(nv), as.integer(ne),
      as.integer(em-1), as.integer(edgeWeightVector(graph)),
      as.integer(init.node-1))
# names(ans) <- c("edgeList", "weights")
# ans$nodes <- nodes(graph)
# ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans+1
})

if (!isGeneric("dfs")) 
   setGeneric("dfs", function(graph)standardGeneric("dfs"))

setMethod("dfs", "graph", function(graph) {
 nv <- length(nodes(graph))
 em <- edgeMatrix(graph)
 ne <- ncol(em)
 ans <- .Call("BGL_dfs_D", as.integer(nv), as.integer(ne),
      as.integer(em-1), as.double(edgeWeightVector(graph)))
 names(ans) <- c("discovered", "finish")
 lapply(ans,function(x)x+1)
})


dijkstra.sp <- function(x,init.ind=1) {
    nN <- nodes(x)
    if (is.character(init.ind)) 
        II <- match(init.ind, nN, 0)
    else
        II <- init.ind
    if (II < 1) 
        stop("bad value for init.ind")
    nv <- length(nN)
    if (II > nv) 
        stop(paste("only", nv, "nodes but init.ind is ", II,
            sep = " "))
    em <- edgeMatrix(x)
    ne <- ncol(em)
    eW <- eWV(x,em)
    ans <- .Call("BGL_dijkstra_shortest_paths_D", as.integer(nv), 
        as.integer(ne), as.integer(em-1), as.double(eW), 
        as.integer(II - 1))
    names(ans) <- c("distances", "penult")
    names(ans[[1]]) <- nN
    ans$penult <- ans$penult + 1
    ans[["start"]] <- II
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

connectedComp <- function (g) 
{
    if (length(agrep("solaris", version$platform))==1) return(
		"inoperative under solaris at present; try windows, linux or BSD")
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    eV <- eWV(g, em)
    x<-.Call("BGL_connected_components_U", as.integer(nv), as.integer(ne), 
        as.integer(em-1), as.double(eV))
    split(nodes(g),x+1)
}

strongComp <- function (g) 
{
    if (edgemode(g) == "undirected") stop("only applicable to directed graphs")
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    x <- .Call("BGL_strong_components_D", as.integer(nv), as.integer(ne), 
        as.integer(em-1), as.double(edgeWeightVector(g)))
    split(nodes(g),x+1)
}

edgeConnectivity <- function (g) 
{
    if (length(agrep("solaris", version$platform))==1) return(
		"inoperative under solaris at present; try windows, linux or BSD")
    if (edgemode(g) == "directed") stop("only applicable to undirected graphs")
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    ans <- .Call("BGL_edge_connectivity_U", as.integer(nv), as.integer(ne), 
        as.integer(em-1), as.double(edgeWeightVector(g)))
    mes <- ans[[2]]
    mes <- lapply(mes,function(x,y) y[x+1], nodes(g)) # +1 for zero-based BGL
    list(connectivity=ans[[1]], minDisconSet=mes)
}
