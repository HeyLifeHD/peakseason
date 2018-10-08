clustered_hm = function(mat_list, rcb_pal = "Blues", sample_names = NULL,
                        title_size = 0.8, top_profile = FALSE, zmins = NULL, zmaxs = NULL,
                        scale = FALSE, file_name = NULL, hm_width = NULL, hm_height = 10, return_mats = FALSE){

  cluster_tbl = mat_list$clustered_matrix[,1:8]
  matrices = mat_list$clustered_matrix[,9:ncol(mat_list$clustered_matrix)]
  colnames(cluster_tbl) = c("chr", "start", "end", "name", "signale", "signal", "cluster_k", "signal")
  cluster_tbl[,cluster_k := cluster_k+1]
  message("Cluster sizes:")
  print(cluster_tbl[,.N,cluster_k])

  mat_size = sum(as.numeric(unlist(tstrsplit(x = mat_list$param["size"], split = ":"))))/(as.numeric(mat_list$param["binSize"]))
  mat_indices = seq(0, ncol(matrices), mat_size)
  mat_indices[1] = 1

  matrices = lapply(1:(length(mat_indices)-1), function(i){
              matrices[,mat_indices[i]:mat_indices[i+1]]
            })

  names(matrices) = mat_list$cdata$sample_names
  #matrices$param = mat_list$param
  #matrices$cluster_row_cut = cluster_tbl$cluster_k

  if(return_mats){
    matrices$param = mat_list$param
    matrices$cdata = mat_list$cdata
    matrices$cluster_tbl = cluster_tbl
    return(matrices)
  }

  ncuts = cluster_tbl$cluster_k
  ncuts = rev(table(ncuts))
  if(!is.null(ncuts)){
    ncuts = cumsum(ncuts/sum(ncuts))
  }

  #return(ncuts)

  ncuts_labs = c()
  # for(i in 1:(length(ncuts)-1)){
  #   print(ncuts[i]-(ncuts[i] - ncuts[i+1]))
  # }

  size = as.character(mat_list$param["size"])
  nsamps = length(matrices)

  if(!is.null(zmins)){
    if(length(zmins) != length(matrices)){
      warning("zmins are recycled")
      zmins = rep(zmins, length(matrices))
    }
  }

  if(!is.null(zmaxs)){
    if(length(zmaxs) != length(matrices)){
      warning("zmaxs are recycled")
      zmaxs = rep(zmaxs, length(matrices))
    }
  }

  if(!is.null(sample_names)){
    if(length(sample_names) != length(matrices)){
      stop("Number of sample names should be equal to number of samples in the matrix")
    }else{
      names(mat_list) = sample_names
    }
  }

  xlabs = c(sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 1), 0,
            sapply(strsplit(x = as.character(size), split = ":", fixed = TRUE), "[[", 2))
  xticks = c(0,
             1/(sum(as.integer(xlabs[1]), as.integer(xlabs[3]))/as.integer(xlabs[1])),
             1)

  if(!is.null(file_name)){
    #Derive optimal height
    if(is.null(hm_width)){
      hm_width = length(mat_list) * 3
    }

    pdf(file = paste0(file_name, ".pdf"), width = hm_width, height = hm_height, paper = "special", bg = "white")
    # tiff(file = , bg = "white", type = "cairo",
    #      height = hm_height, width = hm_width, units = "cm", res = 100)

  }

  if(top_profile){
    lo_mat = matrix(data = 1:(nsamps*2), nrow = 2, ncol = nsamps)
    lo = layout(mat = lo_mat, heights = c(1, 5))
  }else{
    lo_mat = matrix(data = 1:nsamps, nrow = 1)
    lo = layout(mat = lo_mat)
  }

  for(i in 1:nsamps){

    if(top_profile){
      plot_profile_mini(plot_dat = profile_dat, index = i, ylims = yl)
    }

    # if(i == nsamps){
    #   par(mar = c(3,2,1.5,1), las=1, tcl=-.25, font.main=4)
    # }else{
    #   par(mar = c(3,2,1.5,1), las=1, tcl=-.25, font.main=4)
    # }
    par(mar = c(3,2,1.5,1), las=1, tcl=-.25, font.main=4)
    hm.dat = matrices[[i]]
    data.table::setDF(x = hm.dat)

    if(is.null(zmins)){
      zmin = 0
    }else{
      zmin = zmins[i]
    }

    if(scale){
      hm.dat = scale(x = t(hm.dat), scale = TRUE)
    }else{
      hm.dat =  t(apply(hm.dat, 2, rev))
    }

    if(is.null(zmaxs)){
      colMax = apply(hm.dat, 1, max, na.rm = TRUE)
      colMax = which(x = colMax == max(colMax), arr.ind = TRUE)[1]
      zmax = ceiling(max(boxplot.stats(unlist(hm.dat[colMax,]))$stats))
    }else{
      zmax = zmaxs[i]
    }

    hmcols = suppressWarnings(RColorBrewer::brewer.pal(11, rcb_pal))
    hmcols = grDevices::colorRampPalette(hmcols)(255)
    hm.dat[hm.dat >= zmax] = zmax


    if(i == nsamps){
      image(hm.dat, axes = FALSE, col = hmcols, useRaster = TRUE,
            zlim = c(zmin, zmax), xlim = c(-0.2, 1.2))
    }else{
      image(hm.dat, axes = FALSE, col = hmcols, useRaster = TRUE,
            zlim = c(zmin, zmax), xlim = c(-0.2, 1))
    }


    #Add legend
    image(x = c(-0.1, -0.05), y = seq(0, 1, length.out = length(hmcols)-1),
          z = matrix(data = 1:(length(hmcols)-1), nrow = 1), add = TRUE, col = hmcols)
    axis(side = 2, at = seq(0, 1, length.out = 5),
         labels = NA,
         line = -0.6, font.axis = 2, yaxp  = c(1.1, 1.2, 3), lwd = 1.5)
    mtext(text = round(seq(zmin, zmax, length.out = 5), digits = 2), side = 2,
          line = -0.1, at = seq(0, 1, length.out = 5), font = 2)

    title(main = names(matrices)[i], cex.main = title_size, font.main = 2)
    axis(side = 1, at = xticks,
         labels = NA, lty = 1, lwd = 1.5,
         font.axis = 2, cex.axis = 1, line = 0.7)
    mtext(text = c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]), side = 1, line = 1.7, at = xticks, font = 2)


    segments(x0 = 0, y0 = ncuts, x1 = 1, y1 = ncuts, col = "black", lwd = 4)
    if(i == nsamps){
      text(x = 1.1, y = ncuts, labels = rev(1:length(ncuts)), cex = 1.2)
    }

    rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1, border = "black", lwd = 2)
  }

  if(!is.null(file_name)){
    dev.off()
  }

  cluster_tbl
}
