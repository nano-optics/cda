
  cda <- new( "Module" )
  cd <- new( "Module" )
  linear <- new( "Module" )
  array <- new( "Module" )


 .onLoad <- function(libname, pkgname){
     loadRcppModules(direct=FALSE)
 }

