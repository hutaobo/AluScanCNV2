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
  seqlevelsStyle(cnv_gr) <- "UCSC"
  feature <- subsetByOverlaps(cnv_gr, recurr_cnv) # 'recurr_cnv' is the GenomicRanges
  feature <- as.data.frame(feature)
  rownames(feature) <- with(feature, paste0(seqnames, ":", start, "-", end))
  feature <- feature[, "cnv", drop = FALSE]
  feature <- data.frame(t(feature), stringsAsFactors = FALSE)
  pred <- predict(fit, feature, type = "class") # 'fit' is the prediction model
  if(pred == "tumor") {
    result <- TRUE
    cat("\n",
        "#############################################\n",
        "# This blood sample is from Cancer patient. #\n",
        "#############################################\n")
  } else if (pred == "control") {
    result <- FALSE
    cat("\n",
        "#################################################\n",
        "# This blood sample is from Non-Cancer patient. #\n",
        "#################################################\n")
  }
  return(result)
}
