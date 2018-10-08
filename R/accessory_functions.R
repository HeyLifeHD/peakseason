make_bed = function(bed, op_dir){
  #bwtool tool requires only three columns
  op_dir = paste0(op_dir, "/")

  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }

  temp_op_bed = tempfile(pattern = "ezchipviz_", tmpdir = op_dir, fileext = ".bed")

  if(is.data.frame(bed)){
    data.table::setDT(x = bed)
    data.table::fwrite(x = bed[,1:3], file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }else if(as.logical(length(grep(pattern = '\\.gz$', x = bed, fixed = FALSE)))){
    if(length(bed) > 1){
      warning("Multiple BED files provided. Using first one..")
      bed = bed[1]
    }
    bed = gsub(pattern = " ", replacement = "\\ ", x = bed, fixed = TRUE) #Change spaces with \ for unix style paths
    message("Compressed bed file. Deflating..")
    cmd = paste0("zcat < ", bed, " | grep -v ^# | cut -f 1-3  > ", temp_op_bed)
    system(command = cmd, intern = TRUE)
  }else{
    if(length(bed) > 1){
      warning("Multiple BED files provided. Using first one..")
      bed = bed[1]
    }
    bed = gsub(pattern = " ", replacement = "\\ ", x = bed, fixed = TRUE) #Change spaces with \ for unix style paths
    cmd = paste0("cat < ", bed, " |  grep -v ^# | cut -f 1-3 > ", temp_op_bed)
    system(command = cmd, intern = TRUE)
  }

  return(temp_op_bed)
}

id2bed = function(ids){
  chr = unlist(data.table::tstrsplit(x = ids, split = ":", keep = 1))
  loci = unlist(data.table::tstrsplit(x = ids, split = ":", keep = 2))
  loci_start = unlist(tstrsplit(x = loci, split = '-', keep = 1))
  loci_end = unlist(tstrsplit(x = loci, split = '-', keep = 2))
  data.table::data.table(chr = chr, start = loci_start, end = loci_end)
}


extend_bed = function(bed, op_dir, startFrom = "tes", up = 2500, down = 2500){
  bed = gsub(pattern = " ", replacement = "\\ ", x = bed, fixed = TRUE) #Change spaces with \ for unix style paths
  op_dir = paste0(op_dir, "/")

  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }


  if(is.data.frame(bed)){
    data.table::setDT(x = bed)
  }else if(as.logical(length(grep(pattern = '\\.gz$', x = bed, fixed = FALSE)))){
    bed = data.table::fread(cmd = paste0("zcat < ", bed), sep = "\t", header = FALSE)
  }else{
    bed = data.table::fread(input = bed, sep = "\t", header = FALSE)
  }

  temp_op_bed = tempfile(pattern = "ezchipviz_", tmpdir = op_dir, fileext = ".bed")
  colnames(bed)[1:3] = c("chr", "start", "end")

  if(startFrom == "tss"){
    bed[,txStart := as.numeric(start) - up]
    bed[,txEnd := as.numeric(start) + down]
  }else if(startFrom == "tes"){
    bed[,txStart := as.numeric(end) - up]
    bed[,txEnd := as.numeric(end) + down]
  }else if(startFrom == "center"){
    bed[, mp := rowMeans(x = bed[,.(as.numeric(start), as.numeric(end))], na.rm = TRUE)]
    bed[,txStart := as.integer(as.numeric(mp) - up)]
    bed[,txEnd := as.integer(as.numeric(mp) + down)]
    bed[, mp := NULL]
  }else{
    stop("startFrom can only be tss, tes, or center")
  }

  data.table::fwrite(x = bed[,.(chr, txStart, txEnd)],
                     file = temp_op_bed,
                     sep = "\t", row.names = FALSE, col.names = FALSE)

  return(temp_op_bed)
}


make_genome_bed = function(genome, op_dir){
  genome.opts = c("hg19", "hg18", "mm10", "mm9")
  if(!genome %in% genome.opts){
    stop("genome can only be hg19, hg18, mm10, mm8")
  }else{
    genome.file = paste0(genome, "_refseq.bed12.gz")
    genome.file = system.file("extdata", genome.file, package = 'ezChIPviz')
    }
  op_dir = paste0(op_dir, "/")
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }

  temp_op_bed = paste0(op_dir, "/", "temp_", genome, ".bed")
  genome = data.table::fread(cmd = paste0("zcat < ", genome.file))
  data.table::fwrite(x = genome[,.(V1, V2, V3, V4, V5, V6)],
                     file = temp_op_bed,
                     sep = "\t", row.names = FALSE, col.names = FALSE)
  # data.table::fwrite(x = rbind(genome[strand %in% "+", .(chrom, txStart, txEnd)], genome[strand %in% "-", .(chrom, txEnd, txStart)]),
  #                    file = temp_op_bed,
  #                    sep = "\t", row.names = FALSE, col.names = FALSE)

  return(temp_op_bed)
}

# Summarize matrix by mean or median
summarizeMats = function(mats = NULL, summarizeBy = 'mean', group = NULL, collapse_reps = FALSE){

  if(!is.null(group)){
    group_u = unique(group)
    if(collapse_reps){
      summarizedMats = lapply(group_u, function(g){
        x = apply(data.table::rbindlist(l = mats[which(group == g)], fill = TRUE, use.names = TRUE), 2, summarizeBy, na.rm = TRUE)
        x
      })
      names(summarizedMats) = group_u
    }else{
      summarizedMats = lapply(mats[1:(length(mats)-1)], function(x){
        if(!is.null(dim(x))){
          x = apply(x, 2, summarizeBy, na.rm = TRUE)
        }
        x
      })
    }
  }else{
    summarizedMats = lapply(mats[1:(length(mats)-1)], function(x){
      if(!is.null(dim(x))){
        x = apply(x, 2, summarizeBy, na.rm = TRUE)
      }
      x
    })
  }

  summarizedMats$param = mats$param
  summarizedMats
}

# estimate tandard deviation for CI
estimateCI = function(mats = NULL, group = NULL, collapse_reps = FALSE){

  if(!is.null(group)){
    if(collapse_reps){
      group_u = unique(group)
      ciMats = lapply(group_u, function(g){
        x = apply(data.table::rbindlist(l = mats[which(group == g)], fill = TRUE, use.names = TRUE), 2, function(y){
          sd(y, na.rm = TRUE)/sqrt(length(y))
        })
        x
      })
      names(ciMats) = group_u
    }else{
      ciMats = lapply(mats[1:(length(mats)-2)], function(x){
        if(!is.null(dim(x))){
          x = apply(x, 2, function(y){
            sd(y, na.rm = TRUE)/sqrt(length(y))
          })
        }
        x
      })
    }
  }else{
    ciMats = lapply(mats[1:(length(mats)-2)], function(x){
      if(!is.null(dim(x))){
        x = apply(x, 2, function(y){
          sd(y, na.rm = TRUE)/sqrt(length(y))
        })
      }
      x
    })
  }

  ciMats
}

# Add legends outside the margin
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

bwt_mats = function(bw, binSize, bed, size, startFrom, op_dir){

  if(!file.exists(bw)){
    stop(paste0(bw, " does not exists"))
  }

  bn = gsub(pattern = "\\.bw$|\\.bigWig$", replacement = "",
            x = basename(bw), ignore.case = TRUE)
  message(paste0("Processing ", bn, ".."))

  bw = gsub(pattern = " ", replacement = "\\ ", x = bw, fixed = TRUE) #Change spaces with \ for unix style paths

  if(startFrom == "center"){
    cmd = paste0("bwtool matrix -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))
  }else if(startFrom == "start"){
    cmd = paste0("bwtool matrix -starts -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))
  }else if(startFrom == "end"){
    cmd = paste0("bwtool matrix -ends -tiled-averages=", binSize, " ", size, " " , bed , " ", bw, " ", paste0(op_dir, "/", bn, ".matrix"))
  }else{
    stop("startFrom can only be: center, start, or end")
  }
  system(command = cmd, intern = TRUE)
  paste0(op_dir, "/", bn, ".matrix")
}


bwt_mats_k = function(bw, binSize, bed, size, startFrom, op_dir, num_k = NULL){

  bn = tempfile(pattern = "ezchipviz_", tmpdir = op_dir, fileext = ".matrix")
  k_centroids = paste0(bn, ".centroids")

  if(startFrom == "center"){
    cmd = paste0("bwtool matrix -keep-bed -tiled-averages=", binSize, " -cluster=", num_k, " ", size, " " , bed , " ", bw, " ", bn)
  }else if(startFrom == "tss"){
    cmd = paste0("bwtool matrix -keep-bed -starts -tiled-averages=", binSize, " -cluster=", num_k," ", size, " " , bed , " ", bw, " ", bn)
  }else if(startFrom == "tes"){
    cmd = paste0("bwtool matrix -keep-bed -ends -tiled-averages=", binSize, " -cluster=", num_k, " ", size, " " , bed , " ", bw, " ", bn)
  }else{
    stop("startFrom can only be: center, tss, or tes")
  }

  system(command = cmd, intern = TRUE)
  bn
}

#convert size to character
check_size = function(size){
  size = as.character(abs(size))[1:2]
  size = paste0(size, collapse = ":")
  size
}

# Order matrices
order_matrix = function(mats, sortBy = NULLs, k = NULL){

  param = mats$param
  cdata = mats$cdata
  mats = mats[1:(length(mats)-2)]

  mats_avg = lapply(mats, function(x){
                    x = as.matrix(as.data.frame(x = x))
                    apply(x, 1, mean, na.rm = TRUE)
              })

  mats_avg = as.data.frame(mats_avg)
  cluster_row_cut = NULL

  if(sortBy == "mean"){
    mats_avg$oall_avg = rowMeans(mats_avg, na.rm = TRUE)
    row_idx = order(mats_avg$oall_avg, decreasing = TRUE)
  }else if(sortBy == "median"){
    mats_avg$oall_avg = apply(mats_avg, 1, median, na.rm = TRUE)
    row_idx = order(mats_avg$oall_avg, decreasing = TRUE)
  }else if(sortBy == "hclust"){
    set.seed(seed = 1024)
    hc = hclust(d = dist(mats_avg))
    mats_avg$hc_order = hc$order
    if(!is.null(k)){
      mats_avg$cluster = cutree(tree = hc, k = k)
      mats_avg$row_idx = 1:nrow(mats_avg)
      mats_avg_spl = split(mats_avg, f = as.factor(as.character(mats_avg$cluster)))
      cluster_row_cut = unlist(lapply(mats_avg_spl, nrow))
      print(cluster_row_cut)
      mats_avg_spl = lapply(mats_avg_spl, function(x){
                      xhc = hclust(d = dist(x[,1:(ncol(x)-3)]))
                      x$row_idx2 = xhc$order
                      x = x[order(x$row_idx2, decreasing = FALSE),, drop = FALSE]
                      x
                    })
      mats_avg_spl = data.table::rbindlist(l = mats_avg_spl, fill = TRUE, use.names = TRUE)
      #cluster_row_cut = cumsum(xhm[,.N, cluster][,N])
      row_idx = mats_avg_spl$row_idx
    }else{
      row_idx = mats_avg$hc_order
    }
  }

  mats = lapply(mats, function(x){
              x = as.data.frame(x = x)
              x = x[row_idx,,drop = FALSE]
              x
            })

  mats$param = param; mats$cdata = cdata
  mats
}

#Plot profile mini version
plot_profile_mini = function(plot_dat, index = 1, ylims = c(0, 2)){
  #size = as.character(plot_dat$param["size"])

  y = plot_dat[[index]]
  x = (1:(length(y)))/length(y)
  par(mar = c(1,2,1,1))
  plot(x, y, type = "l", lwd = 2, axes = FALSE, xlim = c(-0.2, 1), ylim = ylims, xlab = NA, ylab = NA)
  axis(side = 2, at = ylims, lwd = 2, line = -0.6, labels = NA)
  mtext(side = 2, at = ylims, lwd = 2, line = -0.05, text = ylims, cex = 0.8, font = 2, las = 2)
}

#Take refseq genes and extend TSS/TES and extend them 2500bp
make_genome_bed2 = function(genome, op_dir = "./", startFrom = "tss", up = 2500, down = 2500){

  genome.opts = c("hg19", "hg18", "mm10", "mm9")

  if(!genome %in% genome.opts){
    stop("genome can only be hg19, hg18, mm10, mm8")
  }else{
    genome.file = paste0(genome, "_refseq.bed12.gz")
    genome.file = system.file("extdata", genome.file, package = 'ezChIPviz')
  }

  op_dir = paste0(op_dir, "/")
  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }

  temp_op_bed = paste0(op_dir, "/", "temp_", genome, ".bed")
  genome = data.table::fread(cmd = paste0("zcat < ", genome.file))
  genome = rbind(genome[V6 %in% "+", 1:6], genome[V6 %in% "-", .(V1, V3, V2, V4, V5, V6)], use.names = FALSE)


  if(startFrom == "tss"){
    genome[,txStart := as.numeric(V2) - up]
    genome[,txEnd := as.numeric(V2) + down]
  }else if(startFrom == "tes"){
    genome[,txStart := as.numeric(V3) - up]
    genome[,txEnd := as.numeric(V3) + down]
  }else if (startFrom == "center"){
    genome[, mp := rowMeans(x = genome[,.(as.numeric(V2), as.numeric(V3))], na.rm = TRUE)]
    genome[,txStart := as.integer(as.numeric(mp) - up)]
    genome[,txEnd := as.integer(as.numeric(mp) + down)]
    genome[, mp := NULL]
  }else{
    stop("startFrom can only be tss, tes, or center")
  }

  data.table::fwrite(x = genome[,.(V1, txStart, txEnd, V4, V5, V6)],
                     file = temp_op_bed,
                     sep = "\t", row.names = FALSE, col.names = FALSE)

  return(temp_op_bed)
}


# Add legends outside the margin
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


make_genome_bins = function(genome = NULL, op_dir = "./"){
  chr_size = system.file('extdata', "chr_sizes.txt", package = "ezChIPviz")
  chr_size = data.table::fread(cmd = paste0("zcat < ", chr_size))
}
