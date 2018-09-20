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


#' Print LaTex Table
#'
#' Print LaTeX of given table.
#' @param table A data matrix or data frame
#' @param rownames.include A logical scalar for whether to include row names
#' @return A character scalar of LaTeX table
PrintLaTeXTable = function(table, rownames.include=T) {
  print(xtable::xtable(table), rownames.include=rownames.include)
}
