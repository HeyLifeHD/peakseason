estimatePI = function(control_tss, control_gb, treated_tss, treated_gb, pol2PeakAnnot, bn = NULL){

  pol2.peaks = data.table::fread(input = pol2PeakAnnot)
  #Get active transcripts i.e, tx's with pol2 peaks within 1200 us and 300 ds of TSS
  active.tx = unique(ifelse(test = pol2.peaks$Distance_to_TSS >-1200 & pol2.peaks$Distance_to_TSS < 300, yes = pol2.peaks$Nearest_PromoterID, no = NA))
  message(paste0("Active TXs from ChIP: ", length(active.tx)))

  #Prepare control data (i.e, pol2)
  control_tss = data.table::fread(input = control_tss)
  control_tss[,control_tss_signal := sum/size]
  control_gb = data.table::fread(input = control_gb)
  control_gb[,control_gb_signal := sum/size]
  control = merge(control_tss[,.(name, control_tss_signal)], control_gb[,.(name, control_gb_signal)])

  #Prepare treated data (i.e, pol2+inhibitor)
  treated_tss = data.table::fread(input = treated_tss)
  treated_tss[,treated_tss_signal := sum/size]
  treated_gb = data.table::fread(input = treated_gb)
  treated_gb[,treated_gb_signal := sum/size]
  treated = merge(treated_tss[,.(name, treated_tss_signal)], treated_gb[,.(name, treated_gb_signal)])

  #Process only active tx
  control.atx = control[,tx := sapply(strsplit(x = name, split = ':', fixed = TRUE), "[[", 1)][tx %in% active.tx]

  #Read Refseq ID to gene name annotations
  refseq.id2name = data.table::fread(input = "extdata/hg19_refseq_2017-09-12.refGene")
  refseq.id2name[,id := paste0(name, ":", chrom, "-", txStart)]
  #active.tx.ann = refseq.id2name[name %in% active.tx] #Get annotation of active tx

  ##Order by gene name and signal, Remove duplicated gene names (i.e, keeping unique gene with highest pol2 tx) and Remove short tx (<1000 bp)
  control.atx.filt = merge(x = control.atx, refseq.id2name[,.(name, name2, chrom, txStart, txEnd, size = txEnd - txStart)],
                           by.x = 'tx', by.y = 'name', all.x = TRUE)[order(name2, control_tss_signal, decreasing = TRUE)][!duplicated(name2)][!size < 1000]

  message(paste0("After removing short and duplicated TXs: ", nrow(control.atx.filt)))

  #Removing ovelapping tx
  control.atx.filt.spl = split(control.atx.filt, as.factor(as.character(control.atx.filt$chrom)))

  control.atx.filt.noOverlap = c()
  for(i in 1:length(control.atx.filt.spl)){
    x = control.atx.filt.spl[[i]]
    l = nrow(x)
    x = x[order(as.numeric(as.character(txStart)))]
    x$dist = c(0, as.numeric(x[2:l,txStart]) - as.numeric(x[1:l-1,txEnd]))
    x$dup = ifelse(test = x$dist <= 3000, yes = 1, no = NA)
    rowIndex = which(x$dup == 1)-1
    setDF(x)
    x[rowIndex,"dup"] = 1
    setDT(x)
    x = x[!dup %in% 1]
    control.atx.filt.noOverlap = rbind(control.atx.filt.noOverlap, x)
  }

  message(paste0("After removing overlapping TXs: ", nrow(control.atx.filt.noOverlap)))

  piDat = merge(control.atx.filt.noOverlap[,.(name, name2, control_tss_signal, control_gb_signal)], treated, by = 'name')
  piDat = piDat[,tssDiff := log2(control_tss_signal) - log2(treated_tss_signal)][abs(tssDiff) < 1]

  message(paste0("After removing TXs with differences in Pol2 TSS loading: ", nrow(piDat)))

  piDat[,control_pi := control_tss_signal/control_gb_signal][,treated_pi := treated_tss_signal/treated_gb_signal]

  #Plot ECDF of PI
  plot(ecdf(x = log10(piDat$control_pi)), frame.plot = FALSE, col = 'blue', lwd = 2, xlab = NA,
       main = '', las = 1, ylab = NA, cex.axis = 1.2, xlim = c(0, 2.5), border = NA, yaxt = 'n', xaxt = 'n')
  plot(ecdf(x = log10(piDat$treated_pi)), col = 'red', lwd = 2, add = TRUE)
  axis(side = 2, at = seq(0, 1, 0.2), lwd = 2, font.axis =2, font.lab = 2)
  axis(side = 1, at = seq(0, 2.5, 0.5), lwd = 2, font.axis =2, font.lab = 2)
  mtext(text = "Pausing Index", side = 1, font.lab = 2, line = 2.5)
  mtext(text = "CDF", side = 2, font.lab = 2, line = 2.5)

  pdf(file = paste0("analysis/Pol2_PI/", bn, "_PI.pdf"), width = 5, height = 5, bg = 'white', pointsize = 9)
  plot(ecdf(x = log10(piDat$control_pi)), frame.plot = FALSE, col = 'blue', lwd = 2, xlab = NA,
       main = '', las = 1, ylab = NA, cex.axis = 1.2, xlim = c(0, 2.5), border = NA, yaxt = 'n', xaxt = 'n')
  plot(ecdf(x = log10(piDat$treated_pi)), col = 'red', lwd = 2, add = TRUE)
  axis(side = 2, at = seq(0, 1, 0.2), lwd = 2, font.axis =2, font.lab = 2, cex.axis = 1.3)
  axis(side = 1, at = seq(0, 2.5, 0.5), lwd = 2, font.axis =2, font.lab = 2, cex.axis = 1.3)
  mtext(text = "Pausing Index", side = 1, font.lab = 2, line = 2.5)
  mtext(text = "CDF", side = 2, font.lab = 2, line = 2.5)
  dev.off()

  message("Differences in Pol2 Index between Mock and treated, Wilcox-Rank based non-parametric test:")
  wt = wilcox.test(piDat[,treated_pi], piDat[,control_pi])

  print(broom::tidy(x = wt))

  piDat

}
