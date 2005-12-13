separates = function(a, b, S1, g) {
  if( !is.character(a) || !is.character(b) || !is.character(S1) )
    stop("only vectors of node names allowed")
  if( !is(g, "graph") )
    stop("g must be a graph")

  ng = nodes(g)
  if( any(!(a %in% ng)) || any(!(b %in% ng)) || any(!(S1 %in% ng)) )
     stop("arguments must be nodes in the graph")

  ##sanity check
  if ( any(a %in% S1) || any(b %in% S1) || any(a %in% b) )
    stop("a, b and S1 must be disjoint")
  
  left = ng[ !(ng %in% S1)]
  sg = subGraph(left, g)
  ans = johnson.all.pairs.sp(sg)
  sub1 = ans[a, b]
  if( all(!is.finite(sub1)) ) TRUE else FALSE
} 

