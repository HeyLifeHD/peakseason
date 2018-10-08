#' This function plots top ten GO process from GOSEq results
#' @param gores results generated from \code{\link{bw_go}}
#' @param top Default 10
#' @param font_col text font color
#' @param bar_col color for bars
#' @param font_size font size Default 0.75
plotgo = function(gores, top = 10, font_col = "#4D004B", bar_col= "#BDD7E7", font_size = 0.75){

  par(mar = c(3, 1, 2, 1), tck = -.01)
  xl = ceiling(seq(0, max(-log10(gores[1:top, over_represented_pvalue]), na.rm = TRUE), length.out = 4))
  b = barplot(height = -log10(gores[1:top, over_represented_pvalue]), horiz = TRUE, border = NA,
              col = bar_col, axes = FALSE, xlim = c(min(xl), max(xl)))
  text(x = 0.2, y = b, labels = gores[1:top, term], adj = 0, font = 4, cex = font_size, col = font_col)
  axis(side = 1, at = xl, lwd = 2, line = 0.5, cex.axis = 1)
  mtext(text = "-log10(P-value)", side = 1, line = 1.8, font = 1, cex = 1)
}
