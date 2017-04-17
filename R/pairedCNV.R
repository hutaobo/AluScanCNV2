#' Paired CNV calling
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' pairedCNV()

pairedCNV <- function(sample.5k.doc, control.5k.doc, window.size = c("500k", "400k", "300k", "200k", "100k", "50k"), gender = c("NA", "M", "F"), qOutlier = 0.95, output.path = "./", ...) {
  require(data.table)
  sample.name <- sub(".5k.doc", "", basename(sample.5k.doc))
  control.name <- sub(".5k.doc", "", basename(control.5k.doc))
  window.size <- window.size[1]
  gender <- gender[1]
  factor <- get(paste0("factor.", window.size))  # F
  bin <- get(paste0("bin.", window.size))  # FR
  pos <- get(paste0("pos.", window.size))  # FR2

  outlier <- function(x) {
    x[x > quantile(x, qOutlier)] <- quantile(x, qOutlier)
    return(x)
  }  # remove outliers

  sample.5k.read <- read.table(sample.5k.doc)
  if (!grepl("chr", sample.5k.read[1, 1])) {
    sample.5k.read[, 1] <- paste0("chr", sample.5k.read[, 1])
  }
  sample.5k.read <- sample.5k.read[, c(1:3, 6)]
  sample.5k.read <- merge(bin.5k[, 1:3], sample.5k.read, all = TRUE, sort = FALSE)
  sample.5k.read[is.na(sample.5k.read)] <- 0
  colnames(sample.5k.read) <- c("chr", "start", "end", sample.name)
  sample.5k.read[, 4] <- outlier(sample.5k.read[, 4])
  sample.read <- tapply(sample.5k.read[, 4], as.factor(factor$F), sum)

  control.5k.read <- read.table(control.5k.doc)
  if (!grepl("chr", control.5k.read[1, 1])) {
    control.5k.read[, 1] <- paste0("chr", control.5k.read[, 1])
  }
  control.5k.read <- control.5k.read[, c(1:3, 6)]
  control.5k.read <- merge(bin.5k[, 1:3], control.5k.read, all = TRUE, sort = FALSE)
  control.5k.read[is.na(control.5k.read)] <- 0
  colnames(control.5k.read) <- c("chr", "start", "end", control.name)
  control.5k.read[, 4] <- outlier(control.5k.read[, 4])
  control.read <- tapply(control.5k.read[, 4], as.factor(factor$F), sum)



  library(BSgenome.Hsapiens.UCSC.hg19)
  gcContent <- function(regions, ref = BSgenome.Hsapiens.UCSC.hg19) {
    seq <- getSeq(ref, regions)
    gc <- letterFrequency(seq, "GC")
    acgt <- letterFrequency(seq, "ACGT")
    as.vector(ifelse(acgt == 0, NA, gc/acgt))
  }  # from SomaticSignatures
  bin.gr <- GRanges(seqname = as.character(bin$V1), IRanges(start = bin$V2, end = bin$V3))
  GC <- gcContent(bin.gr)  # GC500k

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

  #### About sex chromosome
  chrX <- pos.5k$chr == 23
  chrY <- pos.5k$chr == 24
  autosome <- pos.5k$chr <= 22

  #### male samples

  if (gender == "M") {
    # Glioma sample: WN is tumor sample and WB is blood control; patient is male

    data <- localCNV4Paired(sample.read, control.read, GC, pos, GCmedian = TRUE)

    write.table(data, paste(output.path, sample.name, "-", control.name, ".local.", window.size, ".paired.seg", sep = ""), row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
  }

  #### female samples

  if (gender == "F") {
    femaleChrs5k <- !chrY
    femaleChrs500k <- FR2$chr <= 23

    data <- localCNV4Paired(sample.read[femaleChrs500k], control.read[femaleChrs500k], GC[femaleChrs500k], pos[femaleChrs500k, ], GCmedian = TRUE)

    write.table(data, paste(output.path, sample.name, "-", control.name, ".local.", window.size, ".paired.seg", sep = ""), row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
  }

}
