#' Paired CNV calling
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' pairedCNV()

pairedCNV <- function(sample.5k.doc, control.5k.doc, window.size = c("500k", "400k", "300k", "200k", "100k", "50k"), gender = c("NA", "M", "F"), qOutlier = 0.95) {
  require(data.table)
  sample_name <- sub(".5k.doc", "", basename(sample.5k.doc))
  control_name <- sub(".5k.doc", "", basename(control.5k.doc))

  window.size <- window.size[1]
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

  #### About sex chromosome
  chrX = pos0$chr==23
  chrY = pos0$chr==24
  autosome = pos0$chr<=22


  # F <- get(paste0('factor.', window.size)) FR <- get(paste0('bin.', window.size)) FR2 <- get(paste0('pos.', window.size))

  library(BSgenome.Hsapiens.UCSC.hg19)
  gcContent <- function(regions, ref = BSgenome.Hsapiens.UCSC.hg19) {
    seq <- getSeq(ref, regions)
    gc <- letterFrequency(seq, "GC")
    acgt <- letterFrequency(seq, "ACGT")
    as.vector(ifelse(acgt == 0, NA, gc/acgt))
  }  # from SomaticSignatures
  bin.5k.gr <- GRanges(seqname = as.character(bin.5k$V1), IRanges(start = bin.5k$V2, end = bin.5k$V3))
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
  write.table(data, paste(sample.name, ".local.", window.size, ".paired.seg", sep = ""), col.name = T, row.name = FALSE, quote = FALSE, sep = "\t")



file.path <- "~/Software/AluScanCNV/Perl scripts/"

######################
####	import data	####
######################

window.size <- "500k"
reads.file.path <- "./"
output.path <- "~/Downloads/"

sample.list <- c("6-3-3", "6-2-3", "6-1-4", "6-2-3", "6-1-4", "6-1-4")
control.list <- c("6b-4", "6b-4", "6b-4", "6-3-3", "6-3-3", "6-2-3")

# WES
sample_name <- read.table(paste0(file.path, "sample_names.txt"), stringsAsFactors = FALSE)
sample_name <- as.vector(t(sample_name))
sample_symbol <- c("M", "M", "T")
control_symbol <- c("T", "N", "N")
sample.list <- apply(expand.grid(sample_name, sample_symbol), 1, paste, collapse = "")
control.list <- apply(expand.grid(sample_name, control_symbol), 1, paste, collapse = "")

sample.list <- combn(paste0("X", c(2:5, 7, 9:12)), 2)[2, ]
control.list <- combn(paste0("X", c(2:5, 7, 9:12)), 2)[1, ]

# import Reads data
read = fread(paste(reads.file.path, "Reads.5k.data", sep = ""), sep="\t", header=T, data.table = FALSE)

# import position data
pos0 = fread(paste(file.path, "5k.pos", sep = ""), sep="\t", header=T) #5k , header

# import sample information of all
sampleInfo = read.table(file=paste(file.path, "All.samples.ID.txt", sep = ""),header=T,sep="\t")

# import test sample ID
# sampleid=read.table(file=paste(file.path, "samples_ID.txt", sep = ""), colClasses = c("character"))
sampleid <- sampleInfo[, c("Sample", "Gender")]
rownames(sampleid) <- sampleid[, 1]

# import factor file and related bin file and pos file
F = read.table(file=paste(file.path, "F", window.size, ".factor", sep = ""), header = T, sep="\t")
FR = read.table(file=paste(file.path, window.size, ".bin", sep = ""), sep="\t", header = FALSE) # chr
FR2 = read.table(file=paste(file.path, window.size, ".pos", sep = ""), sep="\t", header = T) # header

# import reference genome
reference.file = paste(file.path, "hg19.fa", sep = "") #pathway to reference genome (downloaded by user)

############################
####	Global paramters	####
############################

m = dim(read)[1]
n = dim(read)[2]
N = n-3

# calculate GC content
Bin5k = GRanges(seqname = as.character(read$chr), IRanges(start = read$start, end = read$end))
GC5k = exomeCopy::getGCcontent(Bin5k, reference.file)

Bin500k = GRanges(seqname = as.character(FR$V1), IRanges(start = FR$V2, end = FR$V3))
GC500k = exomeCopy::getGCcontent(Bin500k, reference.file)

#
readI = read[, 4:n]
samples = match(sampleid[, 1], names(readI))
#reference = match(referenceid$V1,names(readI))
samplesind = rep(FALSE, dim(readI)[2])
samplesind[samples]=TRUE

# set up outlier to filter out the abnormal part
qOutlier = 0.95

####	primary analysis remove outlier
for(i in 1:N){
	tmp = readI[,i]
	readI[tmp > quantile(tmp, qOutlier),i] = quantile(tmp, qOutlier)
}

# set up a threhold of read depth to lower the noise signal if necessary
# readI[readI<20]=0

#### About sex chromosome
chrX = pos0$chr==23
chrY = pos0$chr==24
autosome = pos0$chr<=22

sampleInd = match(colnames(readI), sampleInfo$Sample)
males = sampleInfo$Gender[sampleInd]=="M"
females = !males

#######################################
#######	Localized CNV calling #########
#######################################

#### paired samples

for(i in 1:length(sample.list))
{
	#### male samples

	if(sampleid[sample.list[i], 2] == "M")
	{
		# Glioma sample: WN is tumor sample and WB is blood control; patient is male
		WN5k_r = localCNV4Paired(readI[,sample.list[i]], readI[,control.list[i]], GC5k, pos0, GCmedian=TRUE)
		WB500k = tapply(readI[,control.list[i]],as.factor(F$F),sum)
		WN500k = tapply(readI[,sample.list[i]],as.factor(F$F),sum)
		WN500k_r = localCNV4Paired(WN500k, WB500k, GC500k, FR2, GCmedian=TRUE)

		write.table(WN500k_r, paste(output.path, sample.list[i], "-", control.list[i], ".local.", window.size, ".paired.seg", sep=""), row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
	}

	#### female samples

	if(sampleid[sample.list[i], 2] == "F")
	{
		femaleChrs5k=!chrY
		femaleChrs500k = FR2$chr<=23
		X216T5k_r = localCNV4Paired(readI[which(femaleChrs5k),sample.list[i]], readI[which(femaleChrs5k),control.list[i]], GC5k[which(femaleChrs5k)], pos0[which(femaleChrs5k),], GCmedian=TRUE)
		X216T500k = tapply(readI[,sample.list[i]],as.factor(F$F),sum)
		X216B500k = tapply(readI[,control.list[i]],as.factor(F$F),sum)
		X216T500k_r = localCNV4Paired(X216T500k[femaleChrs500k], X216B500k[femaleChrs500k], GC500k[femaleChrs500k], FR2[femaleChrs500k,], GCmedian=TRUE)

		write.table(X216T500k_r, paste(output.path, sample.list[i], "-", control.list[i], ".local.", window.size, ".paired.seg", sep=""), row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")
	}
}

}
