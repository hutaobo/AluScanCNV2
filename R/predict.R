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
    rCNV <- recurr_cnv # 'recurr_cnv' is the GenomicRanges
  }
  feature <- subsetByOverlaps(cnv_gr, rCNV)
  feature <- as.data.frame(feature)
  rownames(feature) <- with(feature, paste0(seqnames, ":", start, "-", end))
  feature <- feature[, "cnv", drop = FALSE]
  feature <- data.frame(t(feature), stringsAsFactors = FALSE)
  if (is.null(model)) {
    model <- fit # 'fit' is the prediction model
  }
  pred <- caret::predict.train(object = model, newdata = feature, type = "raw")
  if (!return) cat(as.character(pred)) else return(as.character(pred))
}


#' Cross platform feature selection
#' @param
#' @keywords
#' @export
#' @examples
#' featureSelection()

featureSelection <- function(nonCancerListA, CancerListA, nonCancerListB, CancerListB, Cri = 0.2) {
  library(GenomicRanges)
  df2gr <- function(x) {
    if(is.character(x)) {
      x <- seg2CNV(x)
    }
    x <- subset(x, recurrence >= (ncol(x) - 4) * Cri)
    x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  }
  # choose frequency > 'Cri' in both control and tumor
  nonCancerListA <- df2gr(nonCancerListA)
  CancerListA <- df2gr(CancerListA)
  nonCancerListB <- df2gr(nonCancerListB)
  CancerListB <- df2gr(CancerListB)
  # both platform reached significant enrichment
  nonCancer_recurr <- subsetByOverlaps(nonCancerListA, nonCancerListB)
  Cancer_recurr <- subsetByOverlaps(CancerListA, CancerListB)
  # bind, and remove the overlap of cancer and non-cancer
  a <- nonCancer_recurr
  b <- Cancer_recurr
  recurr_cnv <- c(GenomicRanges::setdiff(a, GenomicRanges::intersect(a, b)), GenomicRanges::setdiff(b, GenomicRanges::intersect(a, b)))
}


#' Cross platform feature selection method_2
#' @param
#' @keywords
#' @export
#' @examples
#' featureSelection2()

featureSelection2 <- function(nonCancerListA, CancerListA, nonCancerListB, CancerListB, Cri = 0.25) {
  library(GenomicRanges)
  df2gr <- function(x) {
    if(is.character(x)) {
      x <- seg2CNV(x)
    }
    x$recurrence <- x$recurrence / (ncol(x) - 4)
    return(x)
  }
  #- choose frequency > 'Cri' in both control and tumor
  nonCancerListA <- df2gr(nonCancerListA)
  CancerListA <- df2gr(CancerListA)
  nonCancerListB <- df2gr(nonCancerListB)
  CancerListB <- df2gr(CancerListB)
  #- both platform reached significant enrichment
  control <- merge(nonCancerListA, nonCancerListB, by = c("chr", "start", "end"), sort = FALSE)
  control$recurrence <- pmin(control$recurrence.x, control$recurrence.y)
  cancer <- merge(CancerListA, CancerListB, by = c("chr", "start", "end"), sort = FALSE)
  cancer$recurrence <- pmin(cancer$recurrence.x, cancer$recurrence.y)
  #- bind, and remove the overlap of cancer and non-cancer
  total <- merge(control, cancer, by = c("chr", "start", "end"), sort = FALSE)
  total <- total[, c('chr', 'start', 'end', 'recurrence.x.x', 'recurrence.y.x', 'recurrence.x', 'recurrence.x.y', 'recurrence.y.y', 'recurrence.y')]
  total$recurrence <- total$recurrence.x - total$recurrence.y
  recurr_cnv <- total[abs(total$recurrence) > Cri, ]
  recurr_cnv <- makeGRangesFromDataFrame(recurr_cnv, keep.extra.columns = TRUE)
}
