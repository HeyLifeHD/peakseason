#' Groups samples based on signal. Useful for QC analysis
#' @param bw_summary output from \code{\link{extract_summary}}
#' @param condition a vector of condition for each bigWig file. Samples belonging to same condition are treated as replicates. Default NULL.
#' @param fract Fraction of data used for PCA or clustering. Default 0.25
#' @param pcs Principal components to draw.
#' @export

plot_summary = function(bw_summary = NULL, log_transform = TRUE, pos_const = 1, condition = NULL){

  anno_data = bw_summary$cdata
  ref_col_idx = which(colnames(anno_data) == condition)

  if(length(ref_col_idx) == 0){
    stop(paste0(condition, " not found!"))
  }

  anno_data = anno_data[,c("sample_names", condition), with = FALSE]
  colnames(anno_data)[ncol(anno_data)] = "cond"
  if(is.null(levels(anno_data$sample_names))){
    anno_data[,cond := as.factor(as.character(cond))]
  }

  hm_dat = bw_summary$summaries
  data.table::setDF(x = hm_dat)
  hm_dat = hm_dat[,c("id", bw_summary$cdata[,sample_names])]

  hm_dat = hm_dat[!duplicated(hm_dat$id),,]
  rownames(hm_dat) = hm_dat$id
  hm_dat = hm_dat[,-1]

  if(log_transform){
    hm_dat = log2(hm_dat + pos_const)
  }

  if(plot_type == "box"){
    hm_dat = data.table::melt(data = hm_dat)
    hm_dat = data.table::setDT(x = hm_dat)
    colnames(hm_dat)[1] = "sample_names"

    if(!is.null(condition)){
      #hm_dat = merge(hm_dat, anno_data, by = "sample_names", all.x = TRUE)
      cond_lvls = levels(anno_data$cond)
      cols = RColorBrewer::brewer.pal(name = "Dark2", n = 8)
      cols = c(cols[1:length(cond_lvls)])
      names(cols) = c(cond_lvls)
      anno_data$cols = cols[anno_data$cond]
      anno_data = anno_data[order(anno_data$cond),,]

      samps = anno_data$cond
      names(samps) =anno_data$sample_names

      hm_dat = hm_dat[,anno_data$sample_names]

      yl = round(c(min(hm_dat, na.rm = TRUE), max(hm_dat, na.rm = TRUE)), digits = 2)
      par(mar = c(5, 4, 4, 2))
      b = boxplot(hm_dat, xaxt="n", boxwex=0.5, outline=TRUE, lty=1, lwd = 1.4, outwex=0,
                  staplewex=0, axes = FALSE, border = anno_data$cols, horizontal = FALSE, ylim = yl, outpch = 19, outcex = 0.5)

      axis(side = 2, at = yl, las = 1, font = 2, lwd = 2.2)
      axis(side = 1, at = 1:length(b$names), labels = b$names, tick = FALSE, las = 2, font = 2, line = -1)
      #text(x = 1:length(b$names), y = 0, labels = b$names, srt = 45, adj = 1)
      add_legend("topleft", legend = names(cols), col= cols, ncol = 4, pch = 15, bty = "n")
    } else{
      yl = round(c(min(hm_dat, na.rm = TRUE), max(hm_dat, na.rm = TRUE)), digits = 2)
      par(mar = c(5, 4, 2, 2))
      b = boxplot(hm_dat, xaxt="n", boxwex=0.5, outline=TRUE, lty=1, lwd = 1.4, outwex=0,
                  staplewex=0, axes = FALSE, border = "gray70", horizontal = FALSE, ylim = yl, outpch = 19, outcex = 0.5)

      axis(side = 2, at = yl, las = 1, font = 2, lwd = 2.2)
      axis(side = 1, at = 1:length(b$names), labels = b$names, tick = FALSE, las = 2, font = 2, line = -1)
      #text(x = 1:length(b$names), y = 0, labels = b$names, srt = 45, adj = 1)
    }
  } else if(plot_type == "hm"){

  }








}
