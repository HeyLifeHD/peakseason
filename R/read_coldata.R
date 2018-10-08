#'Prepares meta data table from bigWig files.
#'Output from this fucntion is passed to all downstream functions.
#' @param coldata a data.frame or a tsv file contaning absolute path to bigWig files and other optional associated sample annotations
#' @param files_idx column index containing path to bigWig files. Required.
#' @param sample_idx column index containing sample names for corresponding bigWig files. Optional. Default NULL - creates one from bigWig files.
#' @export

read_coldata = function(coldata = NULL, files_idx = NULL, sample_idx = NULL){

  if(is.null(coldata)){
    stop("Provide a data.frame or a tsv file with paths to bigWig files and other optional associated sample annotations")
  }

  if(is.data.frame(x = coldata)){
    coldata = data.table::setDT(x = coldata)
  }else if(file.exists(coldata)){
    coldata = data.table::fread(input = coldata, sep = "\t")
  }

  if(is.null(files_idx)){
    print(colnames(coldata))
    stop("Provide column index containing paths to bigwig files. See above from columns available in coldata")
  }else{
    colnames(coldata)[files_idx] = "files"
    coldata$files = as.character(coldata$files)
  }

  message("Checking for files..")
  for(i in 1:nrow(coldata)){
    if(!file.exists(as.character(coldata$files)[i])){
      stop(paste0(as.character(coldata$files)[i], " does not exist!"))
    }
  }

  if(is.null(sample_idx)){
    message("Missing column index for sample names. Creating one from bigWig files for now..")
    coldata$sample_names = gsub(pattern = "\\.bw|\\.bigWig", replacement = "", x = basename(path = as.character(coldata$files)), ignore.case = TRUE)
  }else{
    colnames(coldata)[sample_idx] = "sample_names"
    coldata$sample_names = as.character(coldata$sample_names)
    if(nrow(coldata[duplicated(sample_names)]) > 0){
      stop(paste0("Duplicated sample names ", unique(coldata[duplicated(sample_names)][,sample_names])))
    }
  }

  message("Done!")

  coldata
}
