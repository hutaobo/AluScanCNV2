#' Using 'seg' files to make the CNV file
#'
#' @param
#' @keywords
#' @export
#' @examples
#' seg2CNV()

seg2CNV <- function(seg.list, return = c('df', 'gr')) {
  sample_name <- strsplit(basename(seg.list_KG[i]), '[.]')[[1]][1]
  require(data.table)
  for (i in seq_along(seg.list)) {
    sample_name <- strsplit(basename(seg.list[i]), '[.]')[[1]][1]
    file <- read.table(seg.list[i], header = TRUE, stringsAsFactors = FALSE)
    cnv <- data.frame(chr = file$chromosome, start = file$start, end = file$end, cnv = ifelse(file$p.value >= 0.01, "N", ifelse(file$zScore > 0, "G", "L")), stringsAsFactors = FALSE)
    colnames(cnv)[4] <- sample_name
    if (i == 1) {
      final_cnv <- cnv
    } else {
      final_cnv <- merge(final_cnv, cnv, all = TRUE, sort = FALSE)
    }
  }
  final_cnv[is.na(final_cnv)] <- "N"
  # add 'recurrence' column
  final_cnv$recurrence <- rowSums(final_cnv == "G") + rowSums(final_cnv == "L")
  if (!grepl("chr", final_cnv[1, 1])) {
    final_cnv$chr <- paste0("chr", final_cnv$chr)
  }
  final_cnv$chr <- factor(final_cnv$chr, levels = paste0("chr", c(1:22, "X", "Y")))
  final_cnv <- final_cnv[order(final_cnv$chr, final_cnv$start), ]
  final_cnv$chr <- as.character(final_cnv$chr)
  return <- return[1]
  if (return == 'gr') {
    require(GenomicRanges)
    final_cnv <- makeGRangesFromDataFrame(final_cnv, keep.extra.columns = TRUE)
  }
  return(final_cnv)
}
