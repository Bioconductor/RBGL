leloc <- function(x,n=1,p=1,q=1,xoff=-20, yoff=-10){
 s <- getSpline(x$edges[[n]],p)@cPoints
 xs <- rep(NA,length(s))
 ys <- rep(NA,length(s))
 for (i in 1:length(s))
  {
  xs[i] <- s[[i]]@x
  ys[i] <- s[[i]]@y
  }
 list(x=mean(xs)+xoff, y=mean(ys)+yoff)
}

