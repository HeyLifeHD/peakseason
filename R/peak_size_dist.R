#' Plots peak size distribution
#'
#' @param bed_files path to bed files
#' @param sample_names Default NULL. Correspoding sample names
#' @param orderByMedian Defaukt TRUE
#' @param fs font size
#' @param show Default NULL. Can be mean, ,median or N
#' @export

peak_size = function(bed_files, sample_names = NULL, orderByMedian = TRUE, fs = 0.8, show = NULL){

  bed = lapply(bed_files, function(b){
    x = data.table::fread(b)
    x[,size := V3 - V2]
  })

  if(is.null(sample_names)){
    names(bed) = unlist(data.table::tstrsplit(x = basename(bed_files), split = "\\.", keep = 1))
  }else{
    names(bed) = sample_names
  }

  bcol = c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
           RColorBrewer::brewer.pal(n = 8, name = "Accent"))

  bed = data.table::rbindlist(lapply(bed, function(x) x[,.(size)]), idcol = "sample")

  if(orderByMedian){
    sampleOrder = bed[,median(size),sample][order(V1, decreasing = TRUE)][,sample]
    bed$sample = factor(x = bed$sample, levels = sampleOrder)
  }

  par(mar = c(3, 4, 2, 2))
  b = boxplot(size ~ sample, data = bed, xaxt="n", boxwex=0.5, outline=FALSE, lty=1, lwd = 1.4, outwex=0,
              staplewex=0, axes = FALSE, horizontal = FALSE, border = bcol)

  xl = seq(min(b$stats), max(b$stats), length.out = 4)

  axis(side = 2, at = xl, las = 1, font = 1, lwd = 2.2)
  axis(side = 1, at = 1:length(b$names), labels = b$names, tick = FALSE, las = 2, font = 1, line = -1, cex.axis = fs)

  if(!is.null(show)){
    if(show == "median"){
      bm = bed[,median(size), .(sample)]
      axis(side = 3, at = 1:length(b$names), labels = bm$V1, font = 1, tick = FALSE, line = -1.5, las = 2)
    }else if(show == "mean"){
      bm = bed[,mean(size), .(sample)]
      axis(side = 3, at = 1:length(b$names), labels = bm$V1, font = 1, tick = FALSE, line = -1.5, las = 2)
    }else if(show == "N"){
      axis(side = 3, at = 1:length(b$names), labels = b$N, font = 1, tick = FALSE, line = -1.5, las = 2)
  }else{
    warning("showN can only be median, mean or N")
    }
  }
}
