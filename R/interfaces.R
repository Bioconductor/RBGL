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

prim.minST <- function ( g ) 
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g, duplicates=TRUE) # conform with edgeWeights unlisted
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))

    ans <- .Call("BGL_PRIM_U", as.integer(nv), as.integer(ne), 
                 as.integer(em-1), as.integer(eW), PACKAGE="RBGL")

    rownames(ans[[1]]) <- c("from", "to")
    rownames(ans[[2]]) <- c("weight")
    ans[[1]][1,] <- ans[[1]][1,] + 1
    ans[[1]][2,] <- ans[[1]][2,] + 1
    list("edges"=ans[[1]], "weights"=ans[[2]])
}

if (!isGeneric("bfs")) 
  setGeneric("bfs", function( object, node, checkConn=TRUE) 
             standardGeneric("bfs"))

setMethod("bfs",c("graph", "missing", "missing"),
  function( object, node, checkConn=TRUE)
          bfs(object, nodes(object)[1], TRUE))

setMethod("bfs",c("graph", "character", "missing"),
  function( object, node, checkConn=TRUE)
          bfs(object, node, TRUE))

setMethod("bfs",c("graph", "character", "logical"),
  function( object, node, checkConn)
          bfs(object, node, checkConn))

setMethod("bfs",c("graph", "character", "logical"),
          function( object, node, checkConn=TRUE) {
              nodvec <- nodes(object)
              if (is.na(startind <- match(node,nodvec)))
                  stop("starting node not found in nodes of graph")
              if (checkConn)
              {
                  if (length(connectedComp(object))>1)
                      stop("graph is not connected")
              }
              nv <- length(nodvec)
              em <- edgeMatrix(object,duplicates=TRUE)
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
   setGeneric("dfs", function(object,node,checkConn=TRUE)
              standardGeneric("dfs"))

setMethod("dfs",c("graph", "missing", "missing"),
          function( object, node, checkConn=TRUE)
          dfs(object, nodes(object)[1], TRUE))

setMethod("dfs",c("graph", "character", "missing"),
          function( object, node, checkConn=TRUE)
          dfs(object, node, TRUE))

setMethod("dfs",c("graph", "character", "logical"),
          function( object, node, checkConn=TRUE) {
          nodvec <- nodes(object)
          if (is.na(startind <- match(node,nodvec)))
          {
              warning("starting node not found in nodes of graph,\nnodes element 1 used")
              startind <- 1
          }
          if (checkConn)
            {
            if (length(connectedComp(object))>1) 
                stop("graph is not connected")
            }
          nv <- length(nodvec)
          em <- edgeMatrix(object,duplicates=TRUE)
          ne <- ncol(em)
          if (startind != 1)  # here we rearrange the node references in edgematrix
          {          # to reflect altered start index
            tem <- em
            em[tem == 1] <- startind
            em[tem == startind] <- 1
          }
          ans <- .Call("BGL_dfs_D", as.integer(nv), as.integer(ne),
               as.integer(em-1), as.integer(rep(1,ne)),
               PACKAGE="RBGL")
          fixup <- function(x) { 
             tm <- x;
             x[tm==(1-1)] <- startind-1
             x[tm==(startind-1)] <- (1-1)
             x
             }
          ans <- lapply(ans, fixup)
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

minCut <- function (g)
{
    if (edgemode(g) == "directed") 
       stop("only applicable to undirected graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)

    ans <- .Call("BGL_min_cut_U", as.integer(nv), as.integer(ne),
                 as.integer(em-1), as.double(rep(1.,ne)), PACKAGE="RBGL")

    s_names <- sapply(ans[[2]]+1, function(x) { nodes(g)[x] })
    vs_names <- sapply(ans[[3]]+1, function(x) { nodes(g)[x] })

    if ( length(s_names) > length(vs_names) )
    {  temp <- s_names; s_names <- vs_names; vs_names <- temp }

    list(mincut=ans[[1]], "S"=s_names, "V-S"=vs_names)
}

highlyConnSG <- function (g, sat=3, ldv=c(3, 2, 1))
{
    lldv <- length(ldv)
    x = ldv[1:lldv-1] - ldv[2:lldv]
   
    if ( length(ldv) <= 1 ||
         length(ldv[ldv>0]) != length(ldv) || 
         length(x[x>0]) != length(x) )
       stop("ldv has to be decreasing sequence of positive integers")
    
    if ( sat <= 0 ) stop("sat has to be positive")

    if (edgemode(g) == "directed")
       stop("only applicable to undirected graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)

    ans <- .Call("BGL_highly_conn_sg", as.integer(nv), as.integer(ne),
                 as.integer(em-1), as.double(rep(1.,ne)), 
                 as.integer(sat), as.integer(lldv), as.integer(ldv),
                 PACKAGE="RBGL")

    ans_names <- sapply(ans, function(x) { nodes(g)[x] })
    list(clusters=ans_names)
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

johnson.all.pairs.sp <- function (g) 
{
    nv <- length(nodes(g))
    if (edgemode(g) == "directed") 
        em <- edgeMatrix(g)
    else em <- edgeMatrix(g, TRUE)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))
    ans <- .Call("BGL_johnson_all_pairs_shortest_paths_D", as.integer(nv), 
        as.integer(ne), as.integer(em - 1), as.double(eW), PACKAGE="RBGL")
    tmp <- matrix(ans, nr = length(nodes(g)))
    dimnames(tmp) <- list(nodes(g), nodes(g))
    tmp[ tmp >= .Machine$double.xmax ] <- Inf
    t(tmp)
}

floyd.warshall.all.pairs.sp <- function (g) 
{
    nv <- length(nodes(g))
    if (edgemode(g) == "directed") 
        em <- edgeMatrix(g)
    else em <- edgeMatrix(g, TRUE)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))
    ans <- .Call("BGL_floyd_warshall_all_pairs_shortest_paths_D", as.integer(nv), 
        as.integer(ne), as.integer(em - 1), as.double(eW), PACKAGE="RBGL")
    tmp <- matrix(ans, nr = length(nodes(g)))
    dimnames(tmp) <- list(nodes(g), nodes(g))
    tmp[ tmp >= .Machine$double.xmax ] <- Inf
    t(tmp)
}

bellman.ford.sp <- function(g, start=nodes(g)[1])
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))

    if ( is.character(start) ) s <- match(start, nodes(g), 0)
    else s <- start 

    if ( s <= 0 || s > nv ) stop("start not found in nodes of g")

    ans <- .Call("BGL_bellman_ford_shortest_paths", as.integer(nv), 
           as.integer(ne), as.integer(em - 1), as.double(eW), 
           as.integer(s-1), PACKAGE="RBGL")

    ans[[2]][ ans[[2]] >= .Machine$double.xmax ] <- Inf
    ans[[3]] <- ans[[3]] + 1

    list("all edges minimized"=ans[[1]], "distance"=ans[[2]], 
         "penult"=ans[[3]], "start"=nodes(g)[s])
}

dag.sp <- function(g, start=nodes(g)[1])
{
    if (edgemode(g) == "undirected") 
       stop("only applicable to directed graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))

    if ( is.character(start) ) s <- match(start, nodes(g), 0)
    else s <- start

    if ( s <= 0 || s > nv ) stop("start not found in nodes of g")

    ans <- .Call("BGL_dag_shortest_paths", as.integer(nv), as.integer(ne), 
            as.integer(em - 1), as.double(eW), as.integer(s-1), PACKAGE="RBGL")
    
    ans[[1]][ ans[[1]] >= .Machine$double.xmax ] <- Inf
    ans[[2]] <- ans[[2]] + 1
    list("distance"=ans[[1]], "penult"=ans[[2]], "start"=nodes(g)[s])
}

transitive.closure <- function (g) 
{
    nv <- length(nodes(g))
    if (edgemode(g) == "directed") 
        em <- edgeMatrix(g)
    else em <- edgeMatrix(g, TRUE)
    ne <- ncol(em)
    ans <- .Call("BGL_transitive_closure_D", as.integer(nv), 
        as.integer(ne), as.integer(em - 1), PACKAGE="RBGL")
    
    v_names <- sapply(ans[[1]]+1, function(x) { nodes(g)[x] })
    ans[[1]] <- v_names

    rownames(ans[[2]]) <- c("from", "to")
    f_names <- sapply(ans[[2]][1,]+1, function(x) { nodes(g)[x] })
    t_names <- sapply(ans[[2]][2,]+1, function(x) { nodes(g)[x] })
    ans[[2]][1,] <- f_names
    ans[[2]][2,] <- t_names

    list("nodes"= ans[[1]], "edges"=ans[[2]])
}

max.flow.internal <- function (g, source, sink, method="Edmunds.Karp")
{
    if (edgemode(g) == "undirected") 
       stop("only applicable to directed graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))

    if ( is.character(source) ) s <- match(source, nodes(g), 0)
    else s <- source
    if ( is.character(sink) ) t <- match(sink, nodes(g), 0)
    else t <- sink

    if ( s <= 0 || s > nv || t <= 0 || t > nv )
       stop("both source and sink need to be nodes in the graph")

    # nodes are numbered from 1 in R graph, but from 0 in BGL graph
    if ( method == "Push.Relabel" )
         ans <- .Call("BGL_push_relabel_max_flow", 
                 as.integer(nv), as.integer(ne), 
                 as.integer(em-1), as.integer(eW), 
                 as.integer(s-1), as.integer(t-1), 
                 PACKAGE="RBGL")
    else  # Edmunds.Karp
         ans <- .Call("BGL_edmunds_karp_max_flow", 
                 as.integer(nv), as.integer(ne), 
                 as.integer(em-1), as.integer(eW), 
                 as.integer(s-1), as.integer(t-1), 
                 PACKAGE="RBGL")

    rownames(ans[[2]]) <- c("from", "to")
    rownames(ans[[3]]) <- c("flow")
    f_names <- sapply(ans[[2]][1,]+1, function(x) { nodes(g)[x] })
    t_names <- sapply(ans[[2]][2,]+1, function(x) { nodes(g)[x] })
    ans[[2]][1,] <- f_names
    ans[[2]][2,] <- t_names
    list("maxflow"=ans[[1]], "edges"=ans[[2]], "flows"=ans[[3]])
}

edmunds.karp.max.flow <- function (g, source, sink)
{
    max.flow.internal(g, source, sink, "Edmunds.Karp")
}

push.relabel.max.flow <- function (g, source, sink)
{
    max.flow.internal(g, source, sink, "Push.Relabel")
}

isomorphism <- function(g1, g2)
{
   nv1 <- length(nodes(g1))
   em1 <- edgeMatrix(g1)
   ne1 <- ncol(em1)

   nv2 <- length(nodes(g2))
   em2 <- edgeMatrix(g2)
   ne2 <- ncol(em2)

   ans <- .Call("BGL_isomorphism", 
	        as.integer(nv1), as.integer(ne1), as.integer(em1-1), 
                as.integer(nv2), as.integer(ne2), as.integer(em2-1), 
                PACKAGE="RBGL")

   list("isomorphism"=ans[[1]])
}

cuthill.mckee.ordering <- function(g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_cuthill_mckee_ordering", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   r_names <- sapply(ans[[1]]+1, function(x) { nodes(g)[x] })

   list("reverse cuthill.mckee.ordering"=r_names,
 	"original bandwidth"=ans[[2]], "new bandwidth"=ans[[3]])
}

sequential.vertex.coloring <- function(g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_sequential_vertex_coloring", 
	        as.integer(nv), as.integer(ne), as.integer(em-1),
                PACKAGE="RBGL")

   names(ans[[2]]) = nodes(g)
   list("no. of colors needed"=ans[[1]], "colors of nodes"=ans[[2]])
}

minDegreeOrdering <- function(g, delta=0)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_min_degree_ordering", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), as.integer(delta),
                PACKAGE="RBGL")

    ip_names <- sapply(ans[[1]]+1, function(x) { nodes(g)[x] })
    p_names <- sapply(ans[[2]]+1, function(x) { nodes(g)[x] })

   list("inverse_permutation"=ip_names, "permutation"=p_names)
}

sloan.ordering <- function(g, w1=1, w2=2)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_sloan_ordering", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), w1, w2,
                PACKAGE="RBGL")

   s_names <- sapply(ans[[1]]+1, function(x) { nodes(g)[x] })

   list("sloan.ordering"=s_names, "bandwidth"=ans[[2]], 
	"profile"=ans[[3]], "maxWavefront"=ans[[4]], 
	"aver.wavefront"=ans[[5]], "rms.wavefront"=ans[[6]])
}

bandwidth <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_bandwidth", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   list("bandwidth"=ans[[1]])
}

gprofile <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_profile", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   list("profile"=ans[[1]])
}

ith.wavefront <- function (g, start=nodes(g)[1])
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

    if ( is.character(start) ) s <- match(start, nodes(g), 0)
    else s <- start

    if ( s <= 0 || s > nv )
       stop("starting node needs to be from the graph")

   ans <- .Call("BGL_ith_wavefront", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.integer(s-1),
                PACKAGE="RBGL")

   list("ith.wavefront"=ans[[1]])
}

maxWavefront <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_max_wavefront", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   list("maxWavefront"=ans[[1]])
}

aver.wavefront <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_aver_wavefront", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   list("aver.wavefront"=ans[[1]])
}

rms.wavefront <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_rms_wavefront", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   list("rms.wavefront"=ans[[1]])
}

init.incremental.components <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_init_incremental_components", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   for ( i in 1:ans[[1]] ) ans[-1][[i]] <- ans[-1][[i]] + 1
   ans
}

incremental.components <- function (g)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_incremental_components", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   for ( i in 1:ans[[1]] ) ans[-1][[i]] <- ans[-1][[i]] + 1
   ans
}

same.component <- function (g, node1, node2)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

    if ( is.character(node1) ) v1 <- match(node1, nodes(g), 0)
    else v1 <- node1

    if ( is.character(node2) ) v2 <- match(node2, nodes(g), 0)
    else v2 <- node2

    if ( v1 <= 0 || v1 > nv || v2 <= 0 || v2 > nv )
       stop("nodes need to be from the graph")

   ans <- .Call("BGL_same_component", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
	        as.integer(v1-1), as.integer(v2-1),
                PACKAGE="RBGL")

   ans
}

circle.layout <- function ( g, radius=1 )
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_circle_layout", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(radius),
                PACKAGE="RBGL")

   rownames(ans[[1]]) <- c("x", "y")
   list("circle.layout"=ans[[1]])
}

kamada.kawai.spring.layout <- function ( g, edge_or_side=1, es_length=1 )
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("BGL_kamada_kawai_spring_layout", 
	       as.integer(nv), as.integer(ne), as.integer(em-1), as.integer(eW),
	       as.logical(edge_or_side), as.double(es_length),
               PACKAGE="RBGL")

   rownames(ans[[1]]) <- c("x", "y")
   list("kamada.kawai.spring.layout"=ans[[1]])
}

brandes.betweenness.centrality <- function ( g )
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("BGL_brandes_betweenness_centrality", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.integer(eW),
                PACKAGE="RBGL")

   list("betweenness.centrality.vertices"=ans[[1]],
   "betweenness.centrality.edges"=ans[[2]],
   "relative.betweenness.centrality.vertices"=ans[[3]],
   "dominance"=ans[[4]]
   )
}

betweenness.centrality.clustering <- function(g, threshold=-1, normalize=T )
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("BGL_betweenness_centrality_clustering", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.integer(eW), as.double(threshold), as.logical(normalize),
                PACKAGE="RBGL")

   rownames(ans[[2]]) <- c("from", "to")
   s_names <- sapply(ans[[2]]+1, function(x) { nodes(g)[x] })
   list("no.of.edges" = ans[[1]], 
        "edges"=s_names,
        "edge.betweenness.centrality"=ans[[3]])
}

biConnComp <- function(g)
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    x<-.Call("BGL_biconnected_components_U", as.integer(nv), as.integer(ne),
        as.integer(em-1), as.double(rep(1,ne)), PACKAGE="RBGL")

    list("no. of biconnected components"= x[[1]],
         "biconnected components" = x[[2]])
}

articulationPoints <- function(g)
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    x<-.Call("BGL_articulation_points_U", as.integer(nv), as.integer(ne),
        as.integer(em-1), as.double(rep(1,ne)), PACKAGE="RBGL")

    a_names <- sapply(x[[2]]+1, function(x) { nodes(g)[x] })
    list("no. of articulation points"= x[[1]],
         "articulation points" = a_names)
}

kingOrdering <- function(g)
{
   list("kingOrdering is not implemented yet")
}

randomGraphLayout<- function(g, width=1, height=1)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_random_layout", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(width), as.double(height),
                PACKAGE="RBGL")

   rownames(ans[[1]]) <- c("x", "y")
   list("randomGraphLayout"=ans[[1]])
}

fruchtermanReingoldForceDirectedLayout<- function(g, width=1, height=1)
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_FRFD_layout", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(width), as.double(height),
                PACKAGE="RBGL")

   rownames(ans[[1]]) <- c("x", "y")
   list("fruchtermanReingoldForceDirectedLayout"=ans[[1]])
}

gursoyAtunLayout <- function(g)
{
   list("gursoyAtunLayout is not implemented yet")
}

astarSearch <- function(g)
{
   list("astarSearch is not implemented yet")
}

