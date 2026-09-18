.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath("www", system.file("app/www", package = pkgname))
}
