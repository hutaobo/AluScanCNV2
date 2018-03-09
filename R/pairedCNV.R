#' Paired CNV calling
#'
#' @description https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' control_doc_path <- system.file("extdata/Breast1_b.5k.doc", package = "AluScanCNV")
#' tumor_doc_path <- system.file("extdata/Breast1_1.5k.doc", package = "AluScanCNV")
#' pairedCNV(control.5k.doc = control_doc_path, sample.5k.doc = tumor_doc_path, window.size = "500k", output.path = "./")

pairedCNV <- function(sample.5k.doc, control.5k.doc, window.size = c("500k", "400k", "300k", "250k", "200k", "100k", "50k"), gender = c("M", "F"), qOutlier = 0.95, output.path = "./", replace = TRUE, ...) {
  # Automatically determine whether 'sample.5k.doc' is file_path or dataframe
  if (is.character(sample.5k.doc)) {
    sample.name <- sub(".5k.doc", "", basename(sample.5k.doc))
    control.name <- sub(".5k.doc", "", basename(control.5k.doc))
  } else {
    sample.name <- colnames(sample.5k.doc)[4]
    control.name <- colnames(control.5k.doc)[4]
  }

  window.size <- window.size[1]
  gender <- gender[1]

  # determine if the output file has already exist
  output_file <- paste(output.path, "/", sample.name, "-", control.name, ".local.", window.size, ".paired.seg", sep = "")
  if (replace == FALSE & file.exists(output_file)) {
    return(NULL)
  }

  factor <- get(paste0("factor.", window.size))  # F
  bin <- get(paste0("bin.", window.size))  # FR
  pos <- get(paste0("pos.", window.size))  # FR2

  outlier <- function(x) {
    x[x > quantile(x, qOutlier)] <- quantile(x, qOutlier)
    return(x)
  }  # remove outliers

  if (is.character(sample.5k.doc)) {
    sample.5k.read <- read.table(sample.5k.doc, stringsAsFactors = FALSE)
    sample.5k.read <- sample.5k.read[, c(1:3, 6)]
  } else {
    sample.5k.read <- sample.5k.doc
  }
  if (!grepl("chr", sample.5k.read[1, 1])) {
    sample.5k.read[, 1] <- paste0("chr", sample.5k.read[, 1])
  }
  sample.5k.read <- merge(bin.5k[, 1:3], sample.5k.read, all = TRUE, sort = FALSE, by.x = colnames(bin.5k)[1:3], by.y = colnames(sample.5k.read)[1:3])
  sample.5k.read[is.na(sample.5k.read)] <- 0
  colnames(sample.5k.read) <- c("chr", "start", "end", sample.name)
  sample.5k.read[, 4] <- outlier(sample.5k.read[, 4])
  sample.read <- tapply(sample.5k.read[, 4], as.factor(factor$F), sum)

  if (is.character(control.5k.doc)) {
    control.5k.read <- read.table(control.5k.doc, stringsAsFactors = FALSE)
    control.5k.read <- control.5k.read[, c(1:3, 6)]
  } else {
    control.5k.read <- control.5k.doc
  }
  if (!grepl("chr", control.5k.read[1, 1])) {
    control.5k.read[, 1] <- paste0("chr", control.5k.read[, 1])
  }
  control.5k.read <- merge(bin.5k[, 1:3], control.5k.read, all = TRUE, sort = FALSE, by.x = colnames(bin.5k)[1:3], by.y = colnames(control.5k.read)[1:3])
  control.5k.read[is.na(control.5k.read)] <- 0
  colnames(control.5k.read) <- c("chr", "start", "end", control.name)
  control.5k.read[, 4] <- outlier(control.5k.read[, 4])
  control.read <- tapply(control.5k.read[, 4], as.factor(factor$F), sum)

  GC <- get(paste0("GC.", window.size))

  localCNV4Paired <- function(Sample, Control, GC, Pos, GCmedian = TRUE) {
    # Both of Sample and Control are
    if (GCmedian) {
      index <- apply(cbind(Sample, Control), 1, prod) > 0
      GCGroups <- cut(GC[index], seq(0, 1, 0.05))
      SampleGC <- tapply(Sample[index], GCGroups, median)
      ControlGC <- tapply(Control[index], GCGroups, median)
      SampleL <- SampleGC[GCGroups]
      ControlL <- ControlGC[GCGroups]
      data <- data.frame(chromosome = Pos[index, "chr"], start = Pos[index, "start"], end = Pos[index, "end"], test = Sample[index], ref = Control[index], GClambdaTest = SampleL, GClambdaRef = ControlL, gc = GCGroups)
      results <- cnv.cal(data)
    } else {
      index <- apply(cbind(Sample, Control), 1, prod) > 0
      GCGroups <- cut(GC[index], seq(0, 1, 0.05))
      SampleL <- mean(Sample[index])
      ControlL <- mean(Control[index])
      data <- data.frame(chromosome = Pos[index, "chr"], start = Pos[index, "start"], end = Pos[index, "end"], test = Sample[index], ref = Control[index], GClambdaTest = SampleL, GClambdaRef = ControlL, gc = GCGroups)
      results <- cnv.cal(data)
    }
  }

  #### male samples
  if (gender == "M") {
    data <- localCNV4Paired(sample.read, control.read, GC, pos, GCmedian = TRUE)
    write.table(data, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }

  #### female samples
  if (gender == "F") {
    femaleChrs500k <- pos$chr <= 23
    data <- localCNV4Paired(sample.read[femaleChrs500k], control.read[femaleChrs500k], GC[femaleChrs500k], pos[femaleChrs500k, ], GCmedian = TRUE)
    write.table(data, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
}
