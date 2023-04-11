#' @title Check if a package is installed and install it if not
#' @description This function checks if a package is installed and install it if not.
#' This function is internally used by the package.
#' @param pkg The name of the package to be loaded.
#' @return TRUE if the package is installed, FALSE otherwise.
#' @importFrom BiocManager install
.requirePackage <- function(pkg){
    if(!(pkg %in% .packages(all.available=TRUE)))
    {
        # Try to install the package
        BiocManager::install(pkg, update=FALSE, ask=FALSE)
    }

    require(pkg, character.only=TRUE, quietly=TRUE)
}