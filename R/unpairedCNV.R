library(data.table)

unpairedCNV <- function(sample.5k.doc, sample.name, window.size = "500k", qOutlier = 0.95) {
  factor <- get(paste0("factor.", window.size))  # F
  bin <- get(paste0("bin.", window.size))  # FR
  pos <- get(paste0("pos.", window.size))  # FR2

  outlier <- function(x) {
    x[x > quantile(x, qOutlier)] <- quantile(x, qOutlier)
    return(x)
  }  # remove outliers

  ref.5k.read <- AluScan.ref.5k.reads
  ref.5k.read[, 4:ncol(ref.5k.read)] <- apply(ref.5k.read[, 4:ncol(ref.5k.read)], 2, outlier)
  ref.read <- apply(ref.5k.read[, 4:ncol(ref.5k.read)], 2, function(x) tapply(x, as.factor(factor$F), sum))

  sample.5k.read <- read.table(sample.5k.doc)
  sample.5k.read <- sample.5k.read[, c(1:3, 6)]
  colnames(sample.5k.read) <- c("chr", "start", "end", sample.name)
  sample.5k.read[, 4] <- outlier(sample.5k.read[, 4])
  sample.5k.read <- merge(bin.5k[, 1:3], sample.5k.read, all = TRUE)
  sample.5k.read[is.na(sample.5k.read)] <- 0
  sample.read <- tapply(sample.5k.read[, 4], as.factor(factor$F), sum)



  # F <- get(paste0('factor.', window.size)) FR <- get(paste0('bin.', window.size)) FR2 <- get(paste0('pos.', window.size))

  library(BSgenome.Hsapiens.UCSC.hg19)
  gcContent <- function(regions, ref = BSgenome.Hsapiens.UCSC.hg19) {
    seq <- getSeq(ref, regions)
    gc <- letterFrequency(seq, "GC")
    acgt <- letterFrequency(seq, "ACGT")
    as.vector(ifelse(acgt == 0, NA, gc/acgt))
  }  # from SomaticSignatures
  bin.5k.gr <- GRanges(seqname = as.character(read$chr), IRanges(start = read$start, end = read$end))
  GC.5k <- gcContent(bin.5k.gr)  # GC5k
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
  data <- localCNV4Pool(sample.read, ref.read, GC, pos, GCmedian = TRUE)
  write.table(data, paste(sample.name, ".local.", window.size, ".unpaired.seg", sep = ""), col.name = T, row.name = FALSE, quote = FALSE, sep = "\t")
}
