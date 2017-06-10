#' Predict the tumor occurrence
#' @param
#' @keywords
#' @export
#' @examples
#' predictTumor()

predictTumor <- function(file_path) {
  require(caret)
  file <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  cnv <- data.frame(chr = file$chromosome, start = file$start, end = file$end, cnv = ifelse(file$p.value >= 0.01, 0, ifelse(file$zScore > 0, 1, -1)), stringsAsFactors = FALSE)
  cnv_gr <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  load("recurr_cnv.RData", verbose = FALSE)  # load the GenomicRanges 'recurr_cnv'
  feature <- subsetByOverlaps(cnv_gr, recurr_cnv)
  feature <- as.data.frame(feature)
  rownames(feature) <- with(feature, paste0("chr", seqnames, ":", start, "-", end))
  feature <- feature[, "cnv", drop = FALSE]
  feature <- data.frame(t(feature), stringsAsFactors = FALSE)
  load("fit.RData", verbose = FALSE)  # load the prediction model 'fit'
  pred <- predict(fit, feature, type = "class")
  if(pred == "tumor") {
    cat("\n",
        "############################################\n",
        "# The blood sample is from Cancer patient. #\n",
        "############################################\n")
    return(TRUE)
  } else if (pred == "control") {
    cat("\n",
        "################################################\n",
        "# The blood sample is from Non-Cancer patient. #\n",
        "################################################\n")
    return(FALSE)
  }
}
