tsort <- function(x) 
{
    if ( !isDirected(x) ) stop("requires directed graph")
    nv <- length(nodes(x))
    em <- edgeMatrix(x)
    ne <- ncol(em)
    ans <- .Call("BGL_tsort_D", 
	as.integer(nv), as.integer(ne), as.integer(em-1), 
	PACKAGE="RBGL")
    if ( any(ans != 0) ) ans <- nodes(x)[ans+1]
    else ans <- character(0)
    ans
}

mstree.kruskal <- function(x) 
{
    nv <- length(nodes(x))
    em <- edgeMatrix(x, duplicates=TRUE) # conform with edgeWeights unlisted
    ne <- ncol(em)
    eW <- unlist(edgeWeights(x))
    ans <- .Call("BGL_KMST_D",
      		as.integer(nv), as.integer(ne),
      		as.integer(em-1), as.double(eW),
                PACKAGE="RBGL")

    ans[[1]] <- apply(ans[[1]], 2, function(x, y) y[x+1], nodes(x))
    rownames(ans[[1]]) <- c("from", "to")
    rownames(ans[[2]]) <- c("weight")
    names(ans) <- c("edgeList", "weights")
    ans$nodes <- nodes(x)
    ans

}

prim.minST <- function ( g ) 
{
    .Deprecated("mstree.prim", "RBGL")
    mstree.prim(g)
}

mstree.prim <- function ( g ) 
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g, duplicates=TRUE) # conform with edgeWeights unlisted
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))

    ans <- .Call("BGL_PRIM_U", 
		as.integer(nv), as.integer(ne), 
                as.integer(em-1), as.double(eW), 
		PACKAGE="RBGL")

    ans[[1]] <- apply(ans[[1]], 2, function(x, y) y[x+1], nodes(g))
    rownames(ans[[1]]) <- c("from", "to")
    rownames(ans[[2]]) <- c("weight")
    names(ans) <- c("edgeList", "weights")
    ans$nodes <- nodes(g)
    ans
}

  setGeneric("bfs", function( object, node, checkConn=TRUE) 
             standardGeneric("bfs"))

setMethod("bfs",c("graph", "missing", "missing"),
  function( object, node, checkConn=TRUE)
          bfs(object, nodes(object)[1], TRUE))

setMethod("bfs",c("graph", "missing", "logical"),
  function( object, node, checkConn=TRUE)
          bfs(object, nodes(object)[1], checkConn))

setMethod("bfs",c("graph", "character", "missing"),
  function( object, node, checkConn=TRUE)
          bfs(object, node, TRUE))

setMethod("bfs",c("graph", "character", "logical"),
  function( object, node, checkConn)
          bfs(object, node, checkConn))

setMethod("bfs",c("graph", "character", "logical"),
    function (object, node = nodes(object)[1], checkConn = TRUE) 
    {
    nodvec <- nodes(object)
    if (!checkConn) 
	warning("API is changing: checkConn is disregarded, connectivity always checked.")
    if (is.na(startind <- match(node, nodvec))) 
        stop("starting node not found in nodes of graph")
    if (length(ccc <- connectedComp(object)) > 1) 
    {
        warning("graph is not connected; returning bfs applied to each connected component")
        alln <- lapply(ccc, function(x) nodes(subGraph(x, object)))
        hasStart <- sapply(alln, function(x) node %in% x)
        def <- lapply(ccc[-which(hasStart)], 
		      function(x) bfs(subGraph(x, object)))
	names(def) <- NULL
        wsta <- bfs(subGraph(ccc[[which(hasStart)]], object), node)
        return(c(list(wsta), def))
    }
    nv <- length(nodvec)
    em <- edgeMatrix(object, duplicates = TRUE)
    ne <- ncol(em)
    ans <- .Call("BGL_bfs_D", 
		as.integer(nv), as.integer(ne), as.integer(em - 1), 
		as.integer(rep(1, ne)), as.integer(startind - 1), 
		PACKAGE = "RBGL")
    sapply((ans + 1), function(x, y) y[x], nodes(object))
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
    	  if (!checkConn) 
	     warning("API is changing: checkConn is disregarded, connectivity always checked.")
          nodvec <- nodes(object)
          if (is.na(startind <- match(node,nodvec)))
          {
              warning("starting node not found in nodes of graph,\nnodes element 1 used")
              startind <- 1
          }
        if (length(ccc <- connectedComp(object)) > 1) 
	{
            warning("graph is not connected; returning dfs applied to each connected component")
            def <- lapply(ccc, function(x) dfs(subGraph(x, object)))
       	    names(def) <- NULL
            return(def)
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
          ans <- .Call("BGL_dfs_D", 
			as.integer(nv), as.integer(ne),
               		as.integer(em-1), as.integer(rep(1,ne)),
               		PACKAGE="RBGL")
          fixup <- function(x) 
	  { 
             tm <- x;
             x[tm==(1-1)] <- startind-1
             x[tm==(startind-1)] <- (1-1)
             x
          }
          ans <- lapply(ans, fixup)
          names(ans) <- c("discovered", "finish")
          lapply(ans, function(x, y) y[x+1], nodes(object))
         })

dijkstra.sp <- function(g,start=nodes(g)[1], eW=unlist(edgeWeights(g))) 
{
    if (!is.character(start)) stop("start must be character")
    if (length(start) !=1 ) stop("start must have length 1")

    II <- match(start, nodes(g), 0)

    if ( II == 0 ) stop("start not found in nodes of g")

    if ( isDirected(g) )
            em <- edgeMatrix(g)
    else
            em <- edgeMatrix(g,TRUE)

    if ( any(eW < 0, na.rm=TRUE) ) 
      stop("'dijkstra.sp' requires that all edge weights are nonnegative")

    nN <- nodes(g)
    nv <- length(nN)
    ne <- ncol(em)
    ans <- .Call("BGL_dijkstra_shortest_paths_D", 
		as.integer(nv), as.integer(ne), 
		as.integer(em-1), as.double(eW), as.integer(II - 1), 
		PACKAGE="RBGL")

    names(ans) <- c("distances", "penult")
    ans[["distances"]][ ans[["distances"]] == .Machine$double.xmax ] <- Inf
    names(ans[["distances"]]) <- names(ans[["penult"]]) <- nN
    ans$penult <- ans$penult + 1
    ans[["start"]] <- II
    names(ans[["start"]]) <- nN[II]
    ans
}

connectedComp <- function (g)
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    x<-.Call("BGL_connected_components_U", 
		as.integer(nv), as.integer(ne),
        	as.integer(em-1), as.double(rep(1,ne)), 
		PACKAGE="RBGL")
    split(nodes(g),x+1)
}

strongComp <- function (g)
{
    if (!isDirected(g)) stop("only applicable to directed graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    x <- .Call("BGL_strong_components_D", 
		as.integer(nv), as.integer(ne),
        	as.integer(em-1), as.double(rep(1,ne)), 
		PACKAGE="RBGL")
    split(nodes(g),x+1)
}

edgeConnectivity <- function (g)
{
    if (isDirected(g)) stop("only applicable to undirected graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    ans <- .Call("BGL_edge_connectivity_U", 
		as.integer(nv), as.integer(ne),
                as.integer(em-1), as.double(rep(1.,ne)), 
		PACKAGE="RBGL")
    mes <- ans[[2]]
    mes <- lapply(mes,function(x,y) y[x+1], nodes(g)) # +1 for zero-based BGL
    list(connectivity=ans[[1]], minDisconSet=mes)
}

minCut <- function (g)
{
    if (isDirected(g)) stop("only applicable to undirected graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)

    ans <- .Call("BGL_min_cut_U", 
		as.integer(nv), as.integer(ne),
                as.integer(em-1), as.double(rep(1.,ne)), 
		PACKAGE="RBGL")

    s_names <- sapply(ans[[2]]+1, function(x) { nodes(g)[x] })
    vs_names <- sapply(ans[[3]]+1, function(x) { nodes(g)[x] })

    if ( length(s_names) > length(vs_names) )
    {  temp <- s_names; s_names <- vs_names; vs_names <- temp }

    list(mincut=ans[[1]], "S"=s_names, "V-S"=vs_names)
}

removeSelfLoops <- function(g)
{
    g1 <- g
    for ( n in nodes(g) )
    {
        if ( n %in% adj(g, n)[[1]] ) g1 <- removeEdge(n, n, g1)
    }
    g1
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

    if (isDirected(g)) stop("only applicable to undirected graphs")

    for ( n in nodes(g) )
       if ( n %in% adj(g, n)[[1]] )
          stop("graph contains self-circle(s), use 'removeSelfLoops' first.")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)

    ans <- .Call("BGL_highly_conn_sg", 
		as.integer(nv), as.integer(ne),
                as.integer(em-1), as.double(rep(1.,ne)), 
                as.integer(sat), as.integer(lldv), as.integer(ldv),
                PACKAGE="RBGL")

    ans_names <- lapply(ans, function(x) { nodes(g)[x] })
    list(clusters=ans_names)
}


extractPath <- function(s, f, pens) {
# use list of penultimates (from dijkstra.sp) to establish
# linear path from node s to node f
    path <- f
    maxl <- length(pens)
    i <- 0
    while (path[1] != s) {
        if (i > maxl)        # no path available
        {
            path <- "NA" 
            break
        }
        path <- c(pens[f], path)
        f <- pens[f]
        i <- i+1
    }
    as.numeric(path)
}

sp.between.scalar <- function (g, start, finish)
{
    stop("sp.between.scalar is obsolete, use sp.between instead")
}

sp.between <- function (g, start, finish, detail=TRUE)
{
    if ( length(start) <= 0 || length(finish) <= 0 )
       stop("missing starting or finishing nodes")

    nG = nodes(g)

    if ( !all(start %in% nG) )
       stop("all starting nodes have to be in the graph")

    if ( !all(finish %in% nG) )
       stop("all finishing nodes have to be in the graph")

    ##get the node index, given the name
    nodeind <- function(n) match(n, nG)

    tmp = cbind(start, finish)
    fl = split(tmp[,2], tmp[,1])

    ustart <- unique(start)
    ans <- list()
    ws <- list()
    eW=edgeWeights(g)
    eWW=unlist(eW)

    if ( any(eWW < 0, na.rm=TRUE) ) 
      stop("'sp.between' requires that all edge weights are nonnegative")

    for (i in 1:length(ustart)) 
    {
        curdi <- dijkstra.sp(g, ustart[i], eWW)$penult
        thiss <- ustart[i]
        thisf <- fl[[thiss]]
        for (j in 1:length(thisf) ) 
        {
             ans[[paste(thiss, thisf[j], sep = ":")]] <-
                   nG[extractPath(nodeind(thiss), nodeind(thisf[j]), curdi)]
        }
    }

    getw <- function(nl) 
    {
         res <- NA

         # obtain weights in g for path of nodes in char vec nl
         if ( length(nl) > 1 )
         {
            res <- rep(NA,length(nl)-1)   # only n-1 pairs
	    wstr <- eW[nl]
            for (i in 1:(length(nl)-1))
               res[i]<-wstr[[i]][nl[i+1]] # use numerical names of weights
	    names(res) <- paste(nl[-length(nl)], nl[-1],
			        sep=ifelse(isDirected(g),"->","--"))
         }
         res
    }
    ws <- lapply(ans, function(x) getw(x))
    lens <- lapply(ws, sum)
    ans2 <- vector("list", length=length(ans))
    names(ans2) = names(ans)
    for (i in 1:length(ans))
    {
       if ( detail ) 
        ans2[[i]] <- list(length=lens[[i]],
                          path_detail=as.vector(ans[[i]]), 
                          length_detail=list(ws[[i]]))
       else
        ans2[[i]] <- list(length=lens[[i]])
    }
    ans2
}

johnson.all.pairs.sp <- function (g) 
{
    nv <- length(nodes(g))
    if (isDirected(g)) 
        em <- edgeMatrix(g)
    else em <- edgeMatrix(g, TRUE)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))
    ans <- .Call("BGL_johnson_all_pairs_shortest_paths_D", 
		as.integer(nv), as.integer(ne), 
		as.integer(em - 1), as.double(eW), 
		PACKAGE="RBGL")
    tmp <- matrix(ans, nr = length(nodes(g)))
    dimnames(tmp) <- list(nodes(g), nodes(g))
    tmp[ tmp >= .Machine$double.xmax ] <- Inf
    t(tmp)
}

floyd.warshall.all.pairs.sp <- function (g) 
{
    nv <- length(nodes(g))
    if (isDirected(g)) 
        em <- edgeMatrix(g)
    else em <- edgeMatrix(g, TRUE)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))
    ans <- .Call("BGL_floyd_warshall_all_pairs_shortest_paths_D", 
		as.integer(nv), as.integer(ne), 
		as.integer(em - 1), as.double(eW), 
		PACKAGE="RBGL")
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

    ans <- .Call("BGL_bellman_ford_shortest_paths", 
		as.integer(nv), as.integer(ne), 
		as.integer(em - 1), as.double(eW), as.integer(s-1), 
		PACKAGE="RBGL")

    ans[[2]][ ans[[2]] >= .Machine$double.xmax ] <- Inf
    ans[[3]] <- ans[[3]] + 1

    names(ans[[2]]) <- names(ans[[3]]) <- nodes(g)

    list("all edges minimized"=ans[[1]], "distance"=ans[[2]], 
         "penult"=ans[[3]], "start"=nodes(g)[s])
}

dag.sp <- function(g, start=nodes(g)[1])
{
    if (!isDirected(g)) stop("only applicable to directed graphs")

    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    eW <- unlist(edgeWeights(g))

    if ( is.character(start) ) s <- match(start, nodes(g), 0)
    else s <- start

    if ( s <= 0 || s > nv ) stop("start not found in nodes of g")

    ans <- .Call("BGL_dag_shortest_paths", 
		as.integer(nv), as.integer(ne), 
            	as.integer(em - 1), as.double(eW), as.integer(s-1), 
		PACKAGE="RBGL")
    
    ans[[1]][ ans[[1]] >= .Machine$double.xmax ] <- Inf
    ans[[2]] <- ans[[2]] + 1

    names(ans[[1]]) <- names(ans[[2]]) <- nodes(g)

    list("distance"=ans[[1]], "penult"=ans[[2]], "start"=nodes(g)[s])
}

transitive.closure <- function (g) 
{
    nv <- length(nodes(g))
    if (isDirected(g)) 
        em <- edgeMatrix(g)
    else em <- edgeMatrix(g, TRUE)
    ne <- ncol(em)
    ans <- .Call("BGL_transitive_closure_D", 
		as.integer(nv), as.integer(ne), as.integer(em - 1), 
		PACKAGE="RBGL")
    
    v_names <- sapply(ans[[1]]+1, function(x) { nodes(g)[x] })
    nv <- length(v_names)

    ans[[2]] <- ans[[2]] + 1
    ne <- ncol(ans[[2]])

    edL <- vector("list", length=nv)
    names(edL) <- v_names
    for(i in 1:nv) edL[[i]] <- list(edges=ans[[2]][2,][which(ans[[2]][1,]==i)])

    g <- new("graphNEL", nodes=v_names, edgeL=edL, edgemode=edgemode(g))

    g

}

max.flow.internal <- function (g, source, sink, method="Edmonds.Karp")
{
    if (!isDirected(g)) stop("only applicable to directed graphs")

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
                 as.integer(em-1), as.double(eW), 
                 as.integer(s-1), as.integer(t-1), 
                 PACKAGE="RBGL")
    else  # Edmonds.Karp
         ans <- .Call("BGL_edmonds_karp_max_flow", 
                 as.integer(nv), as.integer(ne), 
                 as.integer(em-1), as.double(eW), 
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

edmonds.karp.max.flow <- function (g, source, sink)
{
    max.flow.internal(g, source, sink, "Edmonds.Karp")
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
	        as.integer(nv), as.integer(ne), 
		as.integer(em-1), as.integer(s-1),
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

   ans[-1] <- lapply(ans[-1], function(x, y) y[x+1], nodes(g))
   names(ans[[1]]) = "no. of initial components"
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

   ans[-1] <- lapply(ans[-1], function(x, y) y[x+1], nodes(g))
   names(ans[[1]]) = "no. of connected components"
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
   warning("API is changed: use circleLayout instead.")
   circleLayout(g, radius)
}

circleLayout <- function ( g, radius=1 )
{
   if (isDirected(g)) stop("only applicable to undirected graphs")

   if ( radius < 0 ) stop("requires: radius > 0 ")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_circle_layout", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(radius),
                PACKAGE="RBGL")

   rownames(ans) <- c("x", "y")
   colnames(ans) <- nodes(g)
   ans
}

kamada.kawai.spring.layout <- function ( g, edge_or_side=1, es_length=1 )
{
   warning("API is changed: use kamadaKawaiSpringLayout instead.")
   kamadaKawaiSpringLayout(g, edge_or_side, es_length)
}

kamadaKawaiSpringLayout <- function ( g, edge_or_side=1, es_length=1 )
{
   if (isDirected(g)) stop("only applicable to undirected graphs")

   if ( edge_or_side <= 0 ) stop("requires: es_length > 0")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("BGL_kamada_kawai_spring_layout", 
	       as.integer(nv), as.integer(ne), 
	       as.integer(em-1), as.double(eW),
	       as.logical(edge_or_side), as.double(es_length),
               PACKAGE="RBGL")

   rownames(ans) <- c("x", "y")
   colnames(ans) <- nodes(g)
   ans
}

brandes.betweenness.centrality <- function ( g )
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   if ( any(eW <= 0, na.rm=TRUE) ) 
      stop("'brandes.betweenness.centrality' requires that all edge weights are positive")

   ans <- .Call("BGL_brandes_betweenness_centrality", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(eW),
                PACKAGE="RBGL")

   colnames(ans[[1]]) <- nodes(g)
   colnames(ans[[3]]) <- nodes(g)

   rownames(ans[[2]]) <- c("centrality")
   rownames(ans[[5]]) <- c("from", "to")
   ans[[5]] <- apply(ans[[5]], 2, function(x, y) y[x+1], nodes(g))

   list("betweenness.centrality.vertices"=ans[[1]],
        "edges"=ans[[5]],
        "betweenness.centrality.edges"=ans[[2]],
        "relative.betweenness.centrality.vertices"=ans[[3]],
        "dominance"=ans[[4]])
}

betweenness.centrality.clustering <- function(g, threshold=-1, normalize=T )
{
   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("BGL_betweenness_centrality_clustering", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(eW), as.double(threshold), as.logical(normalize),
                PACKAGE="RBGL")

   if ( ans[[1]] > 0 )
   {
      ans[[2]] <- apply(ans[[2]], 2, function(x, y) y[x+1], nodes(g))
      rownames(ans[[2]]) <- c("from", "to")
      rownames(ans[[3]]) <- c("centrality")
   }
   list("no.of.edges" = ans[[1]], 
        "edges"=ans[[2]],
        "edge.betweenness.centrality"=ans[[3]])
}

biConnComp <- function(g)
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    ans <-.Call("BGL_biconnected_components_U", 
		as.integer(nv), as.integer(ne),
        	as.integer(em-1), as.double(rep(1,ne)), 
		PACKAGE="RBGL")

    ans[[2]] <- apply(ans[[2]], 2, function(x, y) y[x+1], nodes(g))
    ans[[3]] <- ans[[3]] + 1    # comp no. starts from 1

    r <- vector("list", ans[[1]])
    for ( i in 1:ans[[1]])
        r[[i]] <- unique(as.vector(ans[[2]][,which(ans[[3]]==i)]))

    r
}

articulationPoints <- function(g)
{
    nv <- length(nodes(g))
    em <- edgeMatrix(g)
    ne <- ncol(em)
    ans <-.Call("BGL_articulation_points_U", 
		as.integer(nv), as.integer(ne),
        	as.integer(em-1), as.double(rep(1,ne)), 
		PACKAGE="RBGL")

    ans <- sapply(ans+1, function(x) { nodes(g)[x] })
    ans
}

kingOrdering <- function(g)
{
   list("kingOrdering is not implemented yet")
}

randomGraphLayout<- function(g, minX=0, maxX=1, minY=0, maxY=1)
{
   if (isDirected(g)) stop("only applicable to undirected graphs")

   if ( minX >= maxX || minY >= maxY )
      stop("requires: minX < maxX and minY < maxY ")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_random_layout", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(minX), as.double(maxX),
		as.double(minY), as.double(maxY),
                PACKAGE="RBGL")

   rownames(ans) <- c("x", "y")
   colnames(ans) <- nodes(g)
   ans
}

fruchtermanReingoldForceDirectedLayout<- function(g, width=1, height=1)
{
   if (isDirected(g)) stop("only applicable to undirected graphs")

   if ( length(connComp(g)) > 1 ) 
      warning("This implementation doesn't handle disconnected graphs well.")

   if ( width <= 0 || height <= 0 )
      stop("requires: width > 0 and height > 0 ")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_FRFD_layout", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.double(width), as.double(height),
                PACKAGE="RBGL")

   rownames(ans) <- c("x", "y")
   colnames(ans) <- nodes(g)
   ans
}

gursoyAtunLayout <- function(g)
{
   list("gursoyAtunLayout is not implemented yet")
}

astarSearch <- function(g)
{
   list("astarSearch is not implemented yet")
}

is.triangulated <- function(g)
{
   if( isDirected(g) ) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("isTriangulated", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   if(ans[1] == 0 ) FALSE else TRUE
}

maxClique <- function(g)
{
   if (isDirected(g)) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("maxClique", 
	        as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   ans_names <- lapply(ans, function(x) { nodes(g)[x] })
   list("maxCliques"=ans_names)
}

kCliques <- function(g)
{
   if( isDirected(g)) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("kCliques", 
		as.integer(nv), as.integer(ne), 
		as.integer(em-1), as.double(eW), 
                PACKAGE="RBGL")

   if ( length(ans) > 0 )
   {
      gn1 <- function(x) { nodes(g)[x+1] }
      gn2 <- function(x) { lapply(x, gn1) }
      ans_names <- lapply(ans, gn2)
      names(ans_names) <- paste(1:length(ans_names), "-cliques", sep="")
   }
   else
      ans_names <- ans

   ans_names
}

lambdaSets <- function(g)
{
   if(isDirected(g)) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
   eW <- unlist(edgeWeights(g))

   ans <- .Call("lambdaSets", 
		as.integer(nv), as.integer(ne), 
		as.integer(em-1), as.double(eW), 
                PACKAGE="RBGL")

   makelist <- function(y)
   {
      z <- table(y) > 1
      z <- names(z[z])
      ans <- vector("list", length=length(z))
      for ( i in 1:length(z) )
          ans[i] <- list(names(y[y==z[i]]))
      list(ans)
   }

   colnames(ans[[2]]) <- nodes(g)
   t <- vector("list", length=nrow(ans[[2]]))
   names(t) <- paste("lambda-", 0:(nrow(ans[[2]])-1), " sets", sep="")
   for ( i in 1:nrow(ans[[2]]) )
       t[i] <- makelist(ans[[2]][i,])

   list("max edge connectivity" = ans[[1]], t)
}

clusteringCoef <- function(g, Weighted=FALSE, vW=degree(g))
{
   if(isDirected(g)) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   if ( nv != length(vW) )
      stop("length(vW) is not equal to number of nodes in the graph")

   ans <- .Call("clusteringCoef", 
		as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.integer(Weighted), as.double(vW), 
                PACKAGE="RBGL")

   list("clustering coefficient" = ans)
}

transitivity <- function(g) 
{
   if(isDirected(g)) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("transitivity", 
		as.integer(nv), as.integer(ne), as.integer(em-1), 
                PACKAGE="RBGL")

   list("transitivity" = ans)
}

clusteringCoefAppr <- function(g, k=length(nodes(g)), Weighted=FALSE, vW=degree(g))
{
   if(isDirected(g)) stop("only appropriate for undirected graphs")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   if ( nv != length(vW) )
      stop("length(vW) is not equal to number of nodes in the graph")

   ans <- .Call("clusteringCoefAppr", as.integer(k), 
		as.integer(nv), as.integer(ne), as.integer(em-1), 
		as.integer(Weighted), as.double(vW), 
                PACKAGE="RBGL")

   list("clustering coefficient" = ans)
}

graphGenerator <- function(n, d, o)
{
   if ( n < 3 )
      stop("Number of nodes (n) should be at least 3.")

   if ( d < 2 )
      stop("Degree of nodes (d) should be at least 2.")

   if ( o <= 0 )
      stop("Parameter (o) should be positive.")

   ans <- .Call("graphGenerator", 
		as.integer(n), as.integer(d), as.integer(o), 
		PACKAGE="RBGL")

   # convert node indexes: from 0-based in C to 1-based in R
   ans[[3]] <- ans[[3]] + 1
   rownames(ans[[3]]) <- c("from", "to")
   
   list("no. of nodes" = ans[[1]], "no.of edges" = ans[[2]], "edges"=ans[[3]])
}

dominatorTree <- function(g, start=nodes(g)[1])
{
   if ( length(start) !=1 ) stop("start must have length 1")

   if ( is.character(start) ) s <- match(start, nodes(g), 0)
   else s <- start 

   nv <- length(nodes(g))
   if ( s <= 0 || s > nv ) stop("start not found in nodes of g")

   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("BGL_dominator_tree", 
                as.integer(nv), as.integer(ne), as.integer(em-1),
                as.integer(s-1),
                PACKAGE="RBGL")

   gn1 <- function(x) { nodes(g)[x+1] }
   ans <- sapply(ans, gn1)
   names(ans) <- nodes(g)
   ans
}

kolmogorovMaxFlow <- function(g)
{
   list("This function is not implemented yet")
}

edmondsMaxCardinalityMatching <- function(g)
{
   list("This function is not implemented yet")
}

sloanStartEndVertices<- function(g)
{
   list("This function is not implemented yet")
}

lengauerTarjanDominatorTree<- function(g)
{
   list("This function is not implemented yet")
}


minimumCycleRatio<- function(g)
{
   list("This function is not implemented yet")
}

maximumCycleRatio<- function(g)
{
   list("This function is not implemented yet")
}


boyerMyrvoldPlanarityTest<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("boyerMyrvoldPlanarityTest",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   ans
}

planarFaceTraversal<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("planarFaceTraversal",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   ans

}

planarCanonicalOrdering<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("planarCanonicalOrdering",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   ans

}

chrobakPayneStraightLineDrawing<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("chrobakPayneStraightLineDrawing",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   ans

}

isStraightLineDrawing<- function(g, drawing)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)
 
   if ( !is.matrix(drawing) || nrow(drawing) != 2 ) #|| ncol(drawing) != nv ) 
      stop("needs 2xm matrix for coordinates")

   ans <- .Call("isStraightLineDrawing",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), as.integer(drawing),
		PACKAGE="RBGL")

   ans

}

isKuratowskiSubgraph<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("isKuratowskiSubgraph",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   ans

}

makeConnected<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("makeConnected",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   #print(ans)

   gn1 <- function(x) { nodes(g)[x+1] }
   e_from = sapply(ans[1, ], gn1)
   e_to   = sapply(ans[2, ], gn1)

   gn = new("graphNEL", nodes=nodes(g), edgemode="undirected")
   gn <- addEdge(e_from, e_to, gn)

   gn 
}

makeBiconnectedPlanar<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("makeBiconnectedPlanar",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   #print(ans)

   if ( ans[[1]] )
   {
      gn1 <- function(x) { nodes(g)[x+1] }
      e_from = sapply(ans[[2]][1, ], gn1)
      e_to   = sapply(ans[[2]][2, ], gn1)

      gn = new("graphNEL", nodes=nodes(g), edgemode="undirected")
      gn = addEdge(e_from, e_to, gn)
   
      gn 
   }
   else
   {
      ans[[1]]
   }
}

makeMaximalPlanar<- function(g)
{
   list("This function is not implemented yet")

   nv <- length(nodes(g))
   em <- edgeMatrix(g)
   ne <- ncol(em)

   ans <- .Call("makeMaximalPlanar",
                as.integer(nv), as.integer(ne),
                as.integer(em-1), PACKAGE="RBGL")

   #print(ans)

   if ( ans[[1]] )
   {
      gn1 <- function(x) { nodes(g)[x+1] }
      e_from = sapply(ans[[2]][1, ], gn1)
      e_to   = sapply(ans[[2]][2, ], gn1)
   
      gn = new("graphNEL", nodes=nodes(g), edgemode="undirected")
      gn = addEdge(e_from, e_to, gn)

      gn
   }
   else
   {
      ans[[1]]
   }

}


