#' Differentail Peak Analysis
#' @details Takes output from \code{extract_summary} and performs differential peak analysis with Limma
#'
#' @param bw_summary Output from \code{\link{extract_summary}}
#' @param condition a column name in \code{coldata} containing sample conditions. Default NULL.
#' @param num Numerator condition. Default NULL
#' @param den Denominator condition. Default NULL
#'
#' @import limma
#' @export
#'

diffpeak = function(bw_summary = NULL, condition = NULL, op_dir = "./", num = NULL, den = NULL){

  exprs = bw_summary$summaries
  coldata = bw_summary$cdata

  data.table::setDF(x = coldata)

  if (is.null(condition)){
    stop("You must define a condition for the differential expression\n")
  }

  if(condition %in% colnames(coldata)){
    condition <- coldata[, condition]
  }else {
    print(coldata)
    stop("The condition argument must be a column name in coldata. See above for list of valid strings")
  }
  design <- model.matrix(~0 + as.factor(condition))
  colnames(design) <- as.character(levels(as.factor(condition)))

  if(is.null(num) & is.null(den)){
    contrast <- vector()
    for (a in 1:length(levels(as.factor(condition)))) {
      for (b in 1:length(levels(as.factor(condition)))) {
        if (a != b) {
          if (a < b) {
            contrast[length(contrast) + 1] <- paste(levels(as.factor(condition))[a],
                                                    levels(as.factor(condition))[b], sep = "-")
          }
        }
      }
    }
    message("Argument num and den are missing. Pefroming diffpeak analysis for below contrast:")
    print(contrast[1])
    num = unlist(data.table::tstrsplit(x = contrast[1], split = "-"))[1]
    den = unlist(data.table::tstrsplit(x = contrast[1], split = "-"))[2]
  }else{
    if(is.null(num) | is.null(den)){
      stop("Num and Den must be provided")
    }else{
      contrast = paste0(num, "-", den)
    }
  }

  exprs[,id := paste0(id, "_spl", 1:nrow(exprs))]
  exprs = exprs[,c("id", coldata$sample_names), with = FALSE]
  data.table::setDF(x = exprs, rownames = exprs$id)
  exprs = exprs[,-1]
  exprs = log2(exprs+1)

  fit <- limma::lmFit(object = exprs, design = design)

  cnt <- paste(colnames(design)[1], colnames(design)[2], sep = "-")
  cMat <- limma::makeContrasts(contrasts = contrast, levels = design)
  fit2 <- limma::contrasts.fit(fit, cMat)
  efit <- limma::eBayes(fit2)

  tt = limma::topTable(fit = efit, coef = 1, number = "all")
  data.table::setDT(x = tt, keep.rownames = TRUE)
  data.table::setDT(bw_summary$cdata)

  #data.table::setDT(x = exprs, keep.rownames = TRUE)
  # print(coldata)
  # cond_samp_names = coldata[which(!coldata[,condition] %in% c(num, den), arr.ind = TRUE), "sample_names",]
  # print(cond_samp_names)

  tt = merge(tt, bw_summary$summaries,
             by.x = "rn", by.y = "id")[order(P.Value, decreasing = FALSE)][order(P.Value, decreasing = FALSE)]
  tt[, rn := unlist(data.table::tstrsplit(rn, split = "_spl", keep = 1))]
  bw_summary$summaries[, id := unlist(data.table::tstrsplit(id, split = "_spl", keep = 1))]
  return(results = tt)
}
