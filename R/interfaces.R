tsort <- function(x) {
 if (edgemode(x) != "directed") stop("requires directed graph")
 nv <- length(nodes(x))
 em <- edgeMatrix(x)
 ne <- ncol(em)
 .Call("BGL_tsort_D", as.integer(nv), as.integer(ne),
      as.integer(em-1), PACKAGE="RBGL")
}

mstree.kruskal <- function(x) {
 nv <- length(nodes(x))
 em <- edgeMatrix(x, duplicates=TRUE) # conform with edgeWeights unlisted
 ne <- ncol(em)
ans <- .Call("BGL_KMST_D",
      		as.integer(nv), as.integer(ne),
      		as.integer(em-1), as.double(unlist(edgeWeights(x))),
             PACKAGE="RBGL")
 names(ans) <- c("edgeList", "weights")
 ans$nodes <- nodes(x)
 ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans
}

if (!isGeneric("bfs")) setGeneric("bfs",
        function( object, node, checkConn=FALSE) standardGeneric("bfs"))

setMethod("bfs",c("graph", "missing", "missing"),
  function( object, node, checkConn=FALSE)
          bfs(object, nodes(object)[1], FALSE))

setMethod("bfs",c("graph", "character", "ANY"),
          function( object, node, checkConn=FALSE) {
              nodvec <- nodes(object)
              if (is.na(startind <- match(node,nodvec)))
                  stop("starting node not found in nodes of graph")
              if (checkConn)
              {
                  if (length(connectedComp(object))>1)
                      stop("graph is not connected")
              }
              nv <- length(nodvec)
              em <- edgeMatrix(object)
              ne <- ncol(em)
              ans <- .Call("BGL_bfs_D", as.integer(nv), as.integer(ne),
                           as.integer(em-1), as.integer(rep(1,ne)),
                           as.integer(startind-1), PACKAGE="RBGL")
              ## names(ans) <- c("edgeList", "weights")
              ## ans$nodes <- nodes(object)
              ## ans[["edgeList"]] <- ans[["edgeList"]] + 1
              ans+1
          })

if (!isGeneric("dfs"))
   setGeneric("dfs", function(object,node,checkConn=FALSE)
              standardGeneric("dfs"))

setMethod("dfs",c("graph", "missing", "missing"),
          function( object, node, checkConn=FALSE)
          dfs(object, nodes(object)[1], FALSE))

setMethod("dfs",c("graph", "character", "ANY"),
          function( object, node, checkConn=FALSE) {
              nodvec <- nodes(object)
              if (is.na(startind <- match(node,nodvec)))
          warning("starting node not found in nodes of graph,\nnodes element 1 used")
 if (node != nodvec[1]) warning("starting node supplied but not equal to nodes(g)[1], which will be used instead.")
 if (checkConn)
   {
   if (length(connectedComp(object))>1) stop("graph is not connected")
   }
 nv <- length(nodvec)
 em <- edgeMatrix(object)
 ne <- ncol(em)
 ans <- .Call("BGL_dfs_D", as.integer(nv), as.integer(ne),
      as.integer(em-1), as.integer(rep(1,ne)),
      PACKAGE="RBGL")
 names(ans) <- c("discovered", "finish")
 lapply(ans,function(x)x+1)
})

dijkstra.sp <- function(g,start=nodes(g)[1]) {
    if (!is.character(start)) stop("start must be character")
    if (length(start) !=1 ) stop("start must have length 1")
    nN <- nodes(g)
    if (is.character(start))
        II <- match(start, nN, 0)
    if (II == 0) stop("start not found in nodes of g")
    nv <- length(nN)
    if (edgemode(g) == "directed")
            em <- edgeMatrix(g)
    else
            em <- edgeMatrix(g,TRUE)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))
    ans <- .Call("BGL_dijkstra_shortest_paths_D", as.integer(nv),
        as.integer(ne), as.integer(em-1), as.double(eW),
        as.integer(II - 1), PACKAGE="RBGL")
    names(ans) <- c("distances", "penult")
    ans[["distances"]][ ans[["distances"]] == .Machine$double.xmax ] <- Inf
    names(ans[["distances"]]) <- names(ans[["penult"]]) <- nN
    ans$penult <- ans$penult + 1
    ans[["start"]] <- II
    names(ans[["start"]]) <- nN[II]
    ans
}

sp.between.old <- function(g, start, finish) {
#
#simple vectorization  of previous sp.between
#
if (any(is.numeric(c(start,finish)))) stop("start and finish are required to be node names; numeric indices not allowed")
#
 if (length(start) == 1) {
  if (length(finish) == 1) return( sp.between.scalar(g, start, finish) )
  else { ans <- lapply( finish, function(x,g,start)
            sp.between.scalar(g,start,x), g=g, start=start )
         names(ans)<-paste(start, finish, sep=":")
         return(ans)
       }
  }
 else if (length(finish) == 1)  {
          ans <- lapply(start,function(x,g,finish)
            sp.between.scalar(g,x,finish), g=g, finish=finish)
          names(ans)  <- paste(start, finish, sep=":")
	  return(ans)
         }
 else if (length(finish) != length(start)) stop("cannot have different nonunity lengths of start and finish")
 else {
       sf <- list();for (i in 1:length(start))  sf[[i]] <- c(start[i],finish[i]);
       ans <- lapply( sf,function(x) sp.between.scalar(g,x[1],x[2]))
          names(ans)  <- paste(start, finish, sep=":")
       return(ans)
      }
 }

sp.between.scalar <- function (g, start, finish)
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
    if (length(f) >1) stop("finish must be scalar")
    if (length(s) >1) stop("start must be scalar")
    no <- nodes(g)
    if (any(is.na(lk <- match(c(s, f), no))))
        stop(paste(paste(c(s, f)[is.na(lk)], collapse = " "),
            "not in nodes of g"))
    s <- (1:length(no))[no == s]
    f <- (1:length(no))[no == f]
    ff <- f
    sp <- dijkstra.sp(g, start)
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
    x<-.Call("BGL_connected_components_U", as.integer(nv), as.integer(ne),
        as.integer(em-1), as.double(rep(1,ne)), PACKAGE="RBGL")
    split(nodes(g),x+1)
}

strongComp <- function (g)
{
    if (edgemode(g) == "undirected") stop("only applicable to directed graphs")
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    x <- .Call("BGL_strong_components_D", as.integer(nv), as.integer(ne),
        as.integer(em-1), as.double(rep(1,ne)), PACKAGE="RBGL")
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
        as.integer(em-1), as.double(rep(1.,ne)), PACKAGE="RBGL")
    mes <- ans[[2]]
    mes <- lapply(mes,function(x,y) y[x+1], nodes(g)) # +1 for zero-based BGL
    list(connectivity=ans[[1]], minDisconSet=mes)
}


extractPath <- function(s, f, pens) {
# use list of penultimates (from dijkstra.sp) to establish
# linear path from node s to node f
    path <- f
    maxl <- length(pens)
    i <- 0
    while (path[1] != s) {
        if (i > maxl)
            stop("penultimates inconsistent with linear path from s to f")
        path <- c(pens[f], path)
        f <- pens[f]
        i <- i+1
    }
    as.numeric(path)
}


sp.between <- function (g, start, finish)
{
    nG = nodes(g)
    ##get the node index, given the name
    nodeind <- function(n) match(n, nG)
    if (length(finish)>=length(start))
        fl <- split(finish, start)
    else if (length(finish)==1)
        fl <- split(rep(finish,length(start)),start)
    ust <- unique(start)
    ans <- list()
    for (i in 1:length(ust)) {
        curdi <- dijkstra.sp(g, ust[i])$penult
        thiss <- ust[i]
        thisf <- fl[[thiss]]
        for (j in 1:length(thisf) ) {
            ans[[paste(thiss, thisf[j], sep = ":")]] <-
                nG[extractPath(nodeind(thiss),
                                   nodeind(thisf[j]), curdi)]
        }
    }
    ws <- list()
    eW = edgeWeights(g)
    getw <- function(nl) {
         # obtain weights in g for path of nodes in char vec nl
	 if (length(nl)<2) stop("sp.between:getw should get paths of length 2 or more")
         res <- rep(NA,length(nl)-1)   # only n-1 pairs
	 wstr <- eW[match(nl, nG)]
         for (i in 1:(length(nl)-1))
            res[i] <-  wstr[[i]][as.character(match(nl[i+1],nG))] # need to use numerical names of weights
	names(res) <- paste(nl[-length(nl)],nl[-1],sep=ifelse(edgemode(g)=="undirected","--","->"))
    res
    }
    ws <- lapply(ans, function(x) getw(x))
    ls <- lapply(ws, sum)
    ans2 <- list()
    ns <- names(ans)
    for (i in 1:length(ns))
      {
      ans2[[ns[i]]] <- list()
      ans2[[ns[i]]]$path <- ans[[ns[i]]]
      ans2[[ns[i]]]$length <- ls[[i]]
      ans2[[ns[i]]]$pweights <- ws[[i]]
      }
    ans2
}
