#' Merge two 'data' files
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' dataMerge()

dataMerge <- function(data_1, data_2) {
  if (!grepl("chr", data1[1, 1])) {
    data1$chr <- paste0("chr", data1$chr)
  }
  if (!grepl("chr", data2[1, 1])) {
    data2$chr <- paste0("chr", data2$chr)
  }
  data <- merge(data1, data2, all = TRUE, sort = FALSE)
  data$chr <- factor(data$chr, levels = paste0("chr", c(1:22, "X", "Y")))
  data <- data[order(data$chr, data$start), ]
  data$chr <- as.character(data$chr)
  return(data)
}


#' Using 'doc' files to make the 'data' file
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' doc2data()

doc2data <- function(doc.list, write.file = TRUE) {
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
  if (write.file == TRUE) {
    write.table(data_file, file = "reads.5k.data", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  } else {
    return(data_file)
  }
}
