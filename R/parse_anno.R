get_geneModel = function(ucsc, gene_name){
  ucsc = data.table::fread(input = ucsc)
  gm = ucsc[geneName %in% gene_name]
  gm = gm[,tx_len := txEnd - txStart][order(tx_len, exonCount, decreasing = TRUE)]

  if(nrow(gm) > 1){
    message("Multiple transcripts found. Using longer isoform.")
    print(gm[,.(geneName, name, tx_len, exonCount)])
    gm = gm[1,]
  }
  gm = gm[,.(geneName, name, chrom, strand, txStart, txEnd, tx_len)]
  exon.tbl = data.table::data.table(start = unlist(strsplit(x = gm[,exonStarts], split = ',')),
                         end = unlist(strsplit(x = gm[,exonEnds], split = ',')))

  return(list(gm, exon.tbl))
}

getOverlaps = function(ucsc = NULL, chr, start, end){
  if(!is.null(ucsc)){
    ucsc = data.table::fread(input = ucsc)
  }else{
    ucsc = system.file("extdata", "ucsc_refseq.gz", package = "ezChIPviz")
    ucsc = data.table::fread(cmd = paste0("zcat <  ", ucsc))
    data.table::setkey(x = ucsc, chrom, txStart, txEnd)
  }

  query = data.table::data.table(chrom= chr, txStart = start, txEnd = end)
  gm = data.table::foverlaps(query, ucsc, type="within")

  if(nrow(gm) > 0){
    gm = gm[,tx_len := txEnd - txStart][order(tx_len, exonCount, decreasing = TRUE)]

    if(nrow(gm) > 1){
      message("Multiple transcripts found. Using longer isoform.")
      print(gm[,.(geneName, name, tx_len, exonCount)])
      gm = gm[1,]
    }
    gm = gm[,.(geneName, name, chrom, strand, txStart, txEnd, tx_len, exonStarts, exonEnds)]
    exon.tbl = data.table::data.table(start = unlist(strsplit(x = gm[,exonStarts], split = ',')),
                                      end = unlist(strsplit(x = gm[,exonEnds], split = ',')))
    gene_model = list(gm, exon.tbl)
  }else{
    gene_model = NULL
  }

  gene_model
}

# plot(x = 0, y = 0, xlim = c(gm[,txStart], gm[,txEnd]), ylim = c(0, 1),
#      axes = FALSE, xlab = '', ylab = '')
# segments(x0 = gm[,txStart], y0 = 0.1, x1 = gm[,txEnd], y1 = 0.1, lwd = 2)
# rect(xleft = exon.tbl$start, ybottom = 0, xright = exon.tbl$end,
#      ytop = 0.2, col = "gray70", border = "gray70")
