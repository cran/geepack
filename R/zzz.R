# .onLoad <- function(libname, pkgname) {
#   library.dynam("geepack", pkgname, libname)
# }

.First.lib <- function(lib, pkg) {
    library.dynam("geepack", pkg, lib)
}

# .onUnload <- function(libpath) {
#   library.dynam.unload("geepack", libpath)
# }
