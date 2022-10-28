# function to be executed on package load
.onLoad <- function(libname, pkgname) {
  assign("pkg.env", new.env(), topenv())
  assign("init", FALSE, pkg.env)
}
