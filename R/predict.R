#' Predict the occurrence of tumor
#' @param
#' @keywords
#' @export
#' @examples
#' predictTumor()

predictTumor <- function(file_path, model = NULL, rCNV = NULL, return = FALSE) {
  require(caret)
  file <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  cnv <- data.frame(chr = file$chromosome, start = file$start, end = file$end, cnv = ifelse(file$p.value >= 0.01, 0, ifelse(file$zScore > 0, 1, -1)), stringsAsFactors = FALSE)
  cnv_gr <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  seqlevelsStyle(cnv_gr) <- "UCSC"
  if(is.null(rCNV)) {
    feature <- subsetByOverlaps(cnv_gr, recurr_cnv) # 'recurr_cnv' is the GenomicRanges
  } else {
    feature <- subsetByOverlaps(cnv_gr, rCNV)
  }
  feature <- as.data.frame(feature)
  rownames(feature) <- with(feature, paste0(seqnames, ":", start, "-", end))
  feature <- feature[, "cnv", drop = FALSE]
  feature <- data.frame(t(feature), stringsAsFactors = FALSE)
  if(is.null(model)) {
    pred <- caret::predict(object = fit, newdata = feature, type = "class") # 'fit' is the prediction model
  } else {
    pred <- caret::predict(object = model, newdata = feature, type = "class")
  }
  result <- ifelse(pred == "tumor", TRUE, FALSE)
  if(!return) {
    if(pred == "tumor") {
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

featureSelection <- function(nonCancerListA, CancerListA, nonCancerListB, CancerListB) {
  if(is.character(nonCancerListA)) {
    nonCancerListA <- doc2data(doc.list = nonCancerListA, write.file = FALSE)
  }
  if(is.character(CancerListA)) {
    CancerListA <- doc2data(doc.list = CancerListA, write.file = FALSE)
  }
  if(is.character(nonCancerListB)) {
    nonCancerListB <- doc2data(doc.list = nonCancerListB, write.file = FALSE)
  }
  if(is.character(CancerListB)) {
    CancerListB <- doc2data(doc.list = CancerListB, write.file = FALSE)
  }
  Cri <- 0.2
  control_recurr <- subset(control_cnv_gr, WGS >= length(WGS_control_blood) * Cri & AluScan >= length(AluScan_control_blood) * Cri)
  tumor_recurr <- subset(tumor_cnv_gr, WGS >= length(WGS_tumor_blood) * Cri & AluScan >= length(AluScan_tumor_blood) * Cri)
  a <- control_recurr[, 0]
  b <- tumor_recurr[, 0]
  recurr_cnv <- c(GenomicRanges::setdiff(a, GenomicRanges::intersect(a, b)), GenomicRanges::setdiff(b, GenomicRanges::intersect(a, b)))
}
