#' Estimates polymerase index from Pol2 ChIP seq data.
#' @param coldata coldata Coldata generated from \code{\link{read_coldata}}
#' @param tx_model Genes in BED12 format.
#' @param tssr_up Region upstream of TSS. Defaulr 300
#' @param tssr_up Region downstream of TSS. Defaulr 50
#' @param gb gene body. Region after Transcription end site. Deafult 3000.
#' @param op_dir Directory to store results. Defult "./"
#' @export


estimate_Pol2PI = function(coldata, tx_model = NULL, tssr_up = 50, tssr_down = 300, gb = 3000, op_dir  = "./"){

  tx_model = data.table::fread(input = tx_model)
  tx_model[,id := paste0(V4, ":", V1, "-", V2)]
  tx_model[,size := as.numeric(V3) - as.numeric(V2)]

  tx_model.neg = tx_model[V6 %in% '-']
  tx_model.pos = tx_model[V6 %in% '+']
  #tx_model.ann = rbind(tx_model.pos, tx_model.neg)
  #tx_model.ann[, id := paste0(V1, ':', V2)]

  #TSSR (–50 bp to +300 bp around TSS)
  tssr = tx_model.pos[,.(V1, V2-tssr_up, V2+tssr_down, V4, id)]
  tssr = rbind(tssr, tx_model.neg[,.(V1, V3-tssr_down, V3+50, V4, id)])
  negative_entries = tssr[V2 < 0 | V3 < 0][, id]

  #Gene body (gene body (+300 downstream of the TSS to +3 kb past the TES)
  gene.body = tx_model.pos[,.(V1, V2+tssr_down, V3+gb, V4, id)] #for genes on watson
  gene.body = rbind(gene.body, tx_model.neg[,.(V1, V2-gb, V3-tssr_down, V4, id)]) #for genes on crick

  negative_entries = c(negative_entries, gene.body[V2 < 0 | V3 < 0][, id])

  if(length(negative_entries) > 0){
    gene.body = gene.body[!id %in% negative_entries]
    tssr = tssr[!id %in% negative_entries]
  }

  data.table::fwrite(x = tssr[,.(V1, V2, V3, id)], file = paste0(op_dir, "/", "tx_model_tssr.bed"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  data.table::fwrite(x = gene.body[,.(V1, V2, V3, id)], file = paste0(op_dir, "/", "tx_model_gb.bed"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  #write.table(gene.body[,.(V1, V2, V3, id)], "extdata/refseq_PI_geneBody.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

  message(paste0("Extracting Pol2 signal around TSS region (+", tssr_up, ", -", tssr_down, ").."))
  tss = suppressMessages(peakSeason::extract_summary(coldata = coldata, bed = paste0(op_dir, "/", "tx_model_tssr.bed"), op_dir = op_dir, remove_dups = FALSE))
  tss = tss[[1]]
  tss[,id := tssr[,id]]
  message(paste0("Extracting Pol2 signal around Genebody region (", gb,  ").."))
  gb = suppressMessages(peakSeason::extract_summary(coldata = coldata, bed = paste0(op_dir, "/", "tx_model_gb.bed"), op_dir = op_dir, remove_dups = FALSE))
  gb = gb[[1]]
  gb[,id := gene.body[,id]]

  pol2_dat = lapply(6:ncol(tss), function(i){
                    #TSSR signal
                    tss_sig = tss[,c(5, 4, i), with = FALSE]
                    colnames(tss_sig)[3] = paste0(colnames(tss_sig)[3], "_tssr")
                    tss_sig[,3] = as.numeric(unlist(x = tss_sig[,3, with = FALSE] )) / as.numeric(unlist(x = tss_sig[,2, with = FALSE] ))
                    tss_sig[,size := NULL]
                    #Gene body signal
                    gb_sig = gb[,c(5, 4, i), with = FALSE]
                    colnames(gb_sig)[3] = paste0(colnames(gb_sig)[3], "_gbr")
                    gb_sig[,3] = as.numeric(unlist(x = gb_sig[,3, with = FALSE] )) / as.numeric(unlist(x = gb_sig[,2, with = FALSE] ))
                    gb_sig[,size := NULL]
                    #Merge tssr and genebody, followed by PI estimation (TSSR/Genebody)
                    colnames(tss_sig)[1] = colnames(gb_sig)[1] = "id"
                    sig = merge(gb_sig, tss_sig, by = "id")
                    sig$pausing_idx =  as.numeric(unlist(x = sig[,3, with = FALSE] )) / as.numeric(unlist(x = sig[,2, with = FALSE] ))
                    sig = merge(sig, tx_model[,.(id, size, V4)])
                    colnames(sig)[ncol(sig)] = "tx_name"
                    #sig = sig[,.(id, tx_name, size, )]
                    #sig[, tx_name := unlist(data.table::tstrsplit(id, split = ":", keep = 1))]
                    sig
                  })

  names(pol2_dat) = coldata$sample_names

  pol2_dat
}
