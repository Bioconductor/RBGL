
.First.lib <- function(libname, pkgname, where) {
    require(graph)
    library.dynam("RBGL", pkgname, libname)
}

