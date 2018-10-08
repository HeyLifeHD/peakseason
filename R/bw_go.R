#' Perform Gene ontology
#' @description This function performs gene ontology for output generatd from  \code{\link{diffpeak}}. For now only works with 'hg19' or 'mm10'
#' @param dp_res results generated from \code{\link{diffpeak}}
#' @param p_val Default 0.01
#' @param fdr Default NULL
#' @param diff_class Default 'all', can be 'up' or 'down'
#' @param bn basename for output tsv file
#' @param plotResults Default TRUE

bw_go = function(dp_res, p_val = 0.01, fdr = NULL, diff_class = "all", bn = NULL, plotResults = TRUE){

  genome = dp_res$param["genome"]
  if(is.na(genome)){
    stop("This funtion works only with diffpeak results from hg19 or mm10 genomes")
  }else{
    message(paste0("Genome: ", genome))
  }


  if(is.null(fdr)){
    go_dat = dp_res$results[P.Value < p_val]
  }else{
    go_dat = dp_res$results[adj.P.Val < fdr]
  }

  if(nrow(go_dat) < 1){
    stop('No diff peaks found')
  }

  if(diff_class == "up"){
    go_dat = go_dat[logFC > 0]
  }else if(diff_class == "down"){
    go_dat = go_dat[logFC < 0]
  }else if(diff_class != "all"){
    stop("diff_class can only be all, up, or down")
  }

  assayedGenes = dp_res$results[,V4]
  deGenes = go_dat[,V4]

  if(genome == "hg19"){
    id2sym = paste0(genome, "_refseq2Gene.tsv.gz")
    id2sym = system.file("extdata", id2sym, package = 'ezChIPviz')
    id2sym = data.table::fread(cmd = paste0("zcat < ", id2sym), header = FALSE)
  }else if(genome == "mm10"){
    id2sym = paste0(genome, "_refseq2Gene.tsv.gz")
    id2sym = system.file("extdata", id2sym, package = 'ezChIPviz')
    id2sym = data.table::fread(cmd = paste0("zcat < ", id2sym), header = FALSE)
  }else{
    stop(paste0("Invalid genome: ",genome, "\nThis funtion works only with diffpeak results from hg19 or mm10 genomes"))
  }

  assayedGenes = unique(id2sym[V1 %in% assayedGenes, V2])
  deGenes = unique(id2sym[V1 %in% deGenes, V2])

  gene.vector = as.integer(assayedGenes %in% deGenes)
  names(gene.vector)= assayedGenes

  pwf = goseq::nullp(gene.vector, genome, "geneSymbol", plot.fit = FALSE)
  GO = goseq::goseq(pwf, genome, "geneSymbol", test.cats = c("GO:BP", "GO:MF"), method = "Hypergeometric")
  data.table::setDT(x = GO)

  GO$fdr = p.adjust(GO$over_represented_pvalue)

  if(!is.null(bn)){
    write.table(GO, paste0(bn, '_GO.tsv'), sep='\t', quote = FALSE, row.names = FALSE)
  }

  if(plotResults){
    plotgo(gores = GO)
  }

  GO
}
