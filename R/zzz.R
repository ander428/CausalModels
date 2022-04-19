# function to be executed on package load
.onLoad <- function(libname, pkgname)
{
  options(warn = -1)
  # library(multcomp)
  assign("pkg.env", new.env(), topenv())
  assign("init", F, pkg.env)
}
