
tsortNEW <- function(x) {
 if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 nv <- length(nodes(x))
 ne <- length(unlist(edges(x)))
 .Call("BGL_tsort_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(x)))
}


KMSTNEW <- function(x) {
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

bfsNEW <- function(x,init.ind=1) {
 if (class(x) != "graphNEL") stop("presently only works for class graphNEL from Bioconductor graph library")
 if (init.ind < 1) stop("use 1-based counting for init.ind")
 nv <- length(nodes(x))
 if (init.ind > nv) stop(paste("only",nv,"nodes but init.ind is ",init.ind,sep=" "))
 ne <- length(unlist(edges(x)))
 ans <- .Call("BGL_BFS_D", as.integer(nv), as.integer(ne),
      as.integer(adjListBGL(x)), as.integer(unlist(edgeWeights(x))),
      as.integer(init.ind-1))
# names(ans) <- c("edgeList", "weights")
# ans$nodes <- nodes(x)
# ans[["edgeList"]] <- ans[["edgeList"]] + 1  # bring to unit-based counting
 ans+1
}
