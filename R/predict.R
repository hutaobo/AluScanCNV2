#' Predict the occurrence of tumor
#' @param
#' @keywords
#' @export
#' @examples
#' cancerPrediction()

cancerPrediction <- function(file_path, model = NULL, rCNV = NULL, return = FALSE) {
  require(caret)
  file <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  cnv <- data.frame(chr = file$chromosome, start = file$start, end = file$end, cnv = ifelse(file$p.value >= 0.01, 0, ifelse(file$zScore > 0, 1, -1)), stringsAsFactors = FALSE)
  cnv_gr <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  seqlevelsStyle(cnv_gr) <- "UCSC"
  # 'feature' is the input of 'predict' function
  if (is.null(rCNV)) {
    feature <- subsetByOverlaps(cnv_gr, recurr_cnv) # 'recurr_cnv' is the GenomicRanges
  } else {
    feature <- subsetByOverlaps(cnv_gr, rCNV)
  }
  feature <- as.data.frame(feature)
  rownames(feature) <- with(feature, paste0(seqnames, ":", start, "-", end))
  feature <- feature[, "cnv", drop = FALSE]
  feature <- data.frame(t(feature), stringsAsFactors = FALSE)
  if (is.null(model)) {
    pred <- caret::predict.train(object = fit, newdata = feature, type = "class") # 'fit' is the prediction model
  } else {
    pred <- caret::predict.train(object = model, newdata = feature, type = "class")
  }
  result <- ifelse(pred == "tumor", TRUE, FALSE)
  if (!return) {
    if (pred == "tumor") {
      cat("\n",
          "#############################################\n",
          "# This blood sample is from Cancer patient. #\n",
          "#############################################\n")
    } else if (pred == "control") {
      cat("\n",
          "#################################################\n",
          "# This blood sample is from Non-Cancer patient. #\n",
          "#################################################\n")
    }
  } else {
    return(result)
  }
}


#' Cross platform feature selection
#' @param
#' @keywords
#' @export
#' @examples
#' featureSelection()

featureSelection <- function(nonCancerListA, CancerListA, nonCancerListB, CancerListB, Cri = 0.2) {
  require(GenomicRanges)
  df2gr <- function(x) {
    if(is.character(x)) {
      x <- seg2CNV(x)
    }
    x <- subset(x, recurrence >= (ncol(x) - 4) * Cri)
    x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  }
  # Choose frequency > 'Cri' in both control and tumor
  nonCancerListA <- df2gr(nonCancerListA)
  CancerListA <- df2gr(CancerListA)
  nonCancerListB <- df2gr(nonCancerListB)
  CancerListB <- df2gr(CancerListB)
  # both platform reached significant enrichment
  nonCancer_recurr <- subsetByOverlaps(nonCancerListA, nonCancerListB)
  Cancer_recurr <- subsetByOverlaps(CancerListA, CancerListB)
  # bind, but remove the overlap of cancer and non-cancer
  a <- nonCancer_recurr
  b <- Cancer_recurr
  recurr_cnv <- c(GenomicRanges::setdiff(a, GenomicRanges::intersect(a, b)), GenomicRanges::setdiff(b, GenomicRanges::intersect(a, b)))
}
