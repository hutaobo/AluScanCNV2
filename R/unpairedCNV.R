#' Unpaired CNV calling
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' unpairedCNV()

unpairedCNV <- function(sample.5k.doc, window.size = c("500k", "400k", "300k", "200k", "100k", "50k"), seq.method = c("AluScan", "WGS"), gender = c("NA", "M", "F"), custom.ref = NULL, custom.ref.info = NULL, qOutlier = 0.95, output.path = "./", doc.file = TRUE, ...) {
  sample.name <- sub(".5k.doc", "", basename(sample.5k.doc))
  window.size <- window.size[1]
  gender <- gender[1]
  factor <- get(paste0("factor.", window.size))  # F
  bin <- get(paste0("bin.", window.size))  # FR
  pos <- get(paste0("pos.", window.size))  # FR2

  outlier <- function(x) {
    x[x > quantile(x, qOutlier)] <- quantile(x, qOutlier)
    return(x)
  }  # remove outliers

  seq.method <- seq.method[1]
  if (is.null(custom.ref)) {
    if (seq.method == "AluScan") {
      ref.5k.read <- AluScan.ref.5k.reads
      ref.info <- AluScan.ref.info
    } else if (seq.method == "WGS") {
      ref.5k.read <- WGS.ref.5k.reads
      ref.info <- WGS.ref.info
    }
  } else {
    ref.5k.read <- custom.ref
    ref.info <- custom.ref.info
  }
  ref.5k.read[, 4:ncol(ref.5k.read)] <- apply(ref.5k.read[, 4:ncol(ref.5k.read)], 2, outlier)
  ref.read <- apply(ref.5k.read[, 4:ncol(ref.5k.read)], 2, function(x) tapply(x, as.factor(factor$F), sum))

  if(doc.file == TRUE) {
    sample.5k.read <- read.table(sample.5k.doc, stringsAsFactors = FALSE)
    sample.5k.read <- sample.5k.read[, c(1:3, 6)]
  }
  if (!grepl("chr", sample.5k.read[1, 1])) {
    sample.5k.read[, 1] <- paste0("chr", sample.5k.read[, 1])
  }
  sample.5k.read <- merge(bin.5k[, 1:3], sample.5k.read, all = TRUE, sort = FALSE)
  sample.5k.read[is.na(sample.5k.read)] <- 0
  colnames(sample.5k.read) <- c("chr", "start", "end", sample.name)
  sample.5k.read[, 4] <- outlier(sample.5k.read[, 4])
  sample.read <- tapply(sample.5k.read[, 4], as.factor(factor$F), sum)

  library(BSgenome.Hsapiens.UCSC.hg19)
  gcContent <- function(regions, ref = BSgenome.Hsapiens.UCSC.hg19) {
    seq <- getSeq(ref, regions)
    gc <- letterFrequency(seq, "GC")
    acgt <- letterFrequency(seq, "ACGT")
    as.vector(ifelse(acgt == 0, NA, gc/acgt))
  }  # from SomaticSignatures
  bin.gr <- GRanges(seqname = as.character(bin$V1), IRanges(start = bin$V2, end = bin$V3))
  GC <- gcContent(bin.gr)  # GC500k

  localCNV4Pool <- function(Test, Ref, GC, Pos, GCmedian = TRUE) {
    if (GCmedian) {
      indexRef <- apply(Ref, 1, prod) > 0
      index <- Test > 0 & indexRef
      refSum <- rowSums(Ref)
      GCGroups <- cut(GC[index], seq(0, 1, 0.05))
      TestGC <- tapply(Test[index], GCGroups, median)
      RefGC <- tapply(refSum[index], GCGroups, median)
      TestL <- TestGC[GCGroups]
      RefL <- RefGC[GCGroups]
      data <- data.frame(chromosome = Pos[index, "chr"], start = Pos[index, "start"], end = Pos[index, "end"], test = Test[index], ref = refSum[index], GClambdaTest = TestL, GClambdaRef = RefL, gc = GCGroups)
      results <- cnv.cal(data)
    } else {
      indexRef <- apply(Ref, 1, prod) > 0
      index <- Test > 0 & indexRef
      refSum <- rowSums(Ref)
      TestL <- mean(Test[index])
      RefL <- mean(refSum[index])
      data <- data.frame(chromosome = Pos[index, "chr"], start = Pos[index, "start"], end = Pos[index, "end"], test = Test[index], ref = refSum[index], GClambdaTest = TestL, GClambdaRef = RefL, gc = GCGroups)
      results <- cnv.cal(data)
    }
  }

  sample.read[sample.read == "chr23"] <- "chrX"
  sample.read[sample.read == "chr24"] <- "chrY"
  if (gender == "NA") {
    ind <- pos$chr <= 22
    data <- localCNV4Pool(sample.read[ind], ref.read[ind, ], GC[ind], pos[ind, ], GCmedian = TRUE)
  } else if (gender == "M") {
    data <- localCNV4Pool(sample.read, ref.read[, ref.info[ref.info$Gender == "M", "Name"]], GC, pos, GCmedian = TRUE)
  } else if (gender == "F") {
    data <- localCNV4Pool(sample.read, ref.read[, ref.info[ref.info$Gender == "F", "Name"]], GC, pos, GCmedian = TRUE)
  }
  write.table(data, paste(output.path, "/", sample.name, ".local.", window.size, ".unpaired.seg", sep = ""), col.name = T, row.name = FALSE, quote = FALSE, sep = "\t")
}
