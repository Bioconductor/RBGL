
.First.lib <- function(libname, pkgname, where) {
    library.dynam("RBGL", pkgname, libname)
    require(methods)
    where <- match(paste("package:", pkgname, sep=""), search())
    .initClasses(where)
    cacheMetaData(as.environment(where))
}

