
.First.lib <- function(libname, pkgname ) {
    require("graph")
    library.dynam("RBGL", pkgname, libname)
}

