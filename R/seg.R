#' Using 'seg' files to make the CNV file
#'
#' @param
#' @keywords
#' @export
#' @examples
#' seg2CNV()

seg2CNV <- function(seg.list) {
  docs <- read.table(doc.list, stringsAsFactors = FALSE)
  for (i in 1:nrow(docs)) {
    require(data.table)
    doc_file <- fread(paste0(dirname(doc.list), "/", docs[i, ]), stringsAsFactors = FALSE, data.table = FALSE)
    sample_name <- strsplit(docs[i, ], split = "[.]")[[1]][1]
    doc_file <- doc_file[, c(1:3, 6)]
    colnames(doc_file) <- c("chr", "start", "end", sample_name)
    if (i == 1) {
      data_file <- doc_file
    } else {
      data_file <- merge(data_file, doc_file, all = TRUE, sort = FALSE)
    }
  }
  if (!grepl("chr", data_file[1, 1])) {
    data_file$chr <- paste0("chr", data_file$chr)
  }
  data_file$chr <- factor(data_file$chr, levels = paste0("chr", c(1:22, "X", "Y")))
  data_file <- data_file[order(data_file$chr, data_file$start), ]
  data_file$chr <- as.character(data_file$chr)
  return(data_file)
}