#' Reload Loaded Package
#'
#' Reload a package already loaded
#' @param pkg A character scalar of package name
#' @return None
ReloadPkg = function(pkg) {
  devtools::reload(pkgload::inst(pkg))
}

#' Detach all Packages
#'
#' Detach all packages from global environment
DetachPkgs = function() {
  nothing = tryCatch(suppressWarnings(lapply(paste('package:',
                                        names(sessionInfo()$otherPkgs),sep=""),
         detach, character.only=TRUE)), error=function(e) {})
}
