#' Checks for bwtool installation binary under system path or user defined path
#' @details By default looks for bwtool under system path.
#' @param path path to bwtool binary. Default looks under system directory.
#' @param warn whether to warn user with log messages. Default TRUE
#'
#' @examples
#' check_bwtools()
#' @export

check_bwtools = function(path = NULL, warn = TRUE){

  if(is.null(path)){
    check = as.character(Sys.which(names = 'bwtool'))[1]
  }else{
    check = paste0(path, "/bwtool")
    if(!file.exists(check)){
        stop("Could not locate bwtool at " , check ,"\nDownload it from here: https://github.com/CRG-Barcelona/bwtool/releases")
    }
  }

  if(check != ""){
    if(warn){
      message(paste0("All good! Found bwtool at: ", check))
    }else{
      return(invisible(0))
    }
  }else{
    stop("Could not locate bwtool. Download it from here: https://github.com/CRG-Barcelona/bwtool/releases")
  }
}
