#' Groups samples based on signal. Useful for QC analysis
#' @param bw_summary output from \code{\link{extract_summary}}
#' @param fract Fraction of data used for PCA or clustering. Default 0.25
#' @param pcs Principal components to draw.
#' @param condition a vector of condition for each bigWig file. Samples belonging to same condition are treated as replicates. Default NULL.
#' @export

bw_pca = function(bw_summary = NULL, condition = NULL, top_fract = 0.25, pcs = c(1, 2)){

  exprs = bw_summary$summaries
  data.table::setDT(bw_summary$cdata)
  coldata = bw_summary$cdata

  pca_dat = exprs[,c(coldata[,sample_names]), with = FALSE]
  data.table::setDF(x = pca_dat)

  pca_dat = log2(pca_dat+1)
  pca_dat$sd = apply(pca_dat, 1, sd)
  pca_dat = pca_dat[order(pca_dat$sd, decreasing = TRUE, na.last = TRUE),,]
  pca_dat = pca_dat[1:(as.integer(top_fract*nrow(pca_dat))), 1:(ncol(pca_dat)-1),]
  pca_dat = pca_dat[complete.cases(pca_dat),,]

  pca_dat_hc = hclust(dist(t(pca_dat)))

  #return(pca_dat)

  p = prcomp(x = pca_dat)
  pdat = as.data.frame(p$rotation[,pcs, drop= FALSE])
  rownames(pdat) = names(p$center)
  var.explained = p$sdev^2/sum(p$sdev^2)

  if(!is.null(condition)){
    data.table::setDF(x = coldata)
    condition = coldata[,condition]
    cols = RColorBrewer::brewer.pal(name = "Dark2", n = 8)
    cols.facs = names(table(condition))
    cols.facs = cols.facs[!cols.facs %in% "NA"]
    cols = c(cols[1:length(cols.facs)], "#B3B3B399")
    names(cols) = c(cols.facs, 'NA')
    cols.fill = cols[condition]
  }else{
    cols.fill = "black"
  }

  yl = c(min(pdat[,2]), max(pdat[,2]))
  xl = c(min(pdat[,1]), max(pdat[,1]))

  layout(mat = matrix(data = c(1, 2, 2), nrow = 3), heights = c(4, 4, 1))
  #par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=4, xpd=TRUE, mar = c(3,4,2,2))
  xlab = paste0("PC", pcs[1] , ": " , round(var.explained[pcs[1]], digits = 2))
  ylab = paste0("PC", pcs[2] , ": " , round(var.explained[pcs[2]], digits = 2))
  plot(pdat, pch = 21, axes = FALSE, xlab = '', ylab = '',
       xlim = xl, ylim = yl, bg = cols.fill, cex = 2, col = "black") #

  #points(pdat, col = grDevices::adjustcolor(col = "black", alpha.f = 0.6))
  mtext(text = ylab, side = 2, line = 2.5, las = 3, font = 2)
  mtext(text = xlab, side = 1, line = 2, las = 1, font = 2)
  axis(side = 1, at = xl, labels = round(xl, digits = 2), lwd = 1.2, font = 2, line = 1.2)
  axis(side = 2, at = yl, labels = round(yl, digits = 2), lwd = 1.2, font = 2, line = 1)

  if(!is.null(condition)){
    legend(x = "top", legend = names(cols), col = cols, pt.lwd = 2,
               bty = "o", border = NA, xpd = TRUE, text.font = 2, pch = 19, ncol = 3)
  }

  plot(pca_dat_hc, hang = -0.5, cex = 0.6, xlab = "", ylab = NA, main = "", lwd = 1.2, cex = 1.2, font = 2)

  data.table::setDT(bw_summary$cdata)

}
