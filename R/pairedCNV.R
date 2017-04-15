#' Paired CNV calling
#'
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4273479/pdf/13336_2014_Article_15.pdf
#' @param
#' @keywords
#' @export
#' @examples
#' pairedCNV()

##################
####	demo.R	####
##################

##################
####	library	####
##################
source("cnvF.R")
library(data.table)
file.path <- "~/Software/AluScanCNV/Perl scripts/"

######################
####	import data	####
######################

window.size <- "500k"
reads.file.path <- "./"
output.path <- "~/Downloads/"

# panel 1 & 2
# sample.list <- c("GC1A", "GC1A", "GC1C", "GC2A", "GC2A", "GC2C", "GC3A", "GC3A", "GC3C", "GC4A", "GC4A", "GC4C", "GC5A", "GC5A", "GC5C",
#	"Breast", "Breast", "Normal", "PBT2T", "PBT2T", "HL3T", "PBT6T", "PBT6T", "HL1T", "PBT7T", "PBT7T", "HL2T", "LC5", "LC5", "LC5N",
#	"LC6", "LC6", "LC6N", "LC11", "LC11", "LC11N", "LC106", "LC106", "LC106N", "LiverT", "LiverT", "LiverN", "LG1T", "LG1T", "LG1N",
#	"LG2T", "LG2T", "LG2N", "LG3T", "LG3T", "LG3N", "LG5T", "LG5T", "LG5N")
#
# control.list <- c("GC1C", "GC1D", "GC1D", "GC2C", "GC2D", "GC2D", "GC3C", "GC3D", "GC3D", "GC4C", "GC4D", "GC4D", "GC5C", "GC5D", "GC5D",
#	"Normal", "Blood", "Blood", "HL3T", "PBT2N", "PBT2N", "HL1T", "PBT6N", "PBT6N", "HL2T", "PBT7N", "PBT7N", "LC5N", "LC5B", "LC5B",
#	"LC6N", "LC6B", "LC6B", "LC11N", "LC11B", "LC11B", "LC106N", "LC106B", "LC106B", "LiverN", "LiverB", "LiverB", "LG1N", "LG1B", "LG1B",
#	"LG2N", "LG2B", "LG2B", "LG3N", "LG3B", "LG3B", "LG5N", "LG5B", "LG5B")

# panel 3
# sample.list <- c("Meta", "Meta", "Meta", "GC1A", "GC1B", "GC1B", "GC2A", "GC2B", "GC2B", "GC3A", "GC3B", "GC3B", "GC4A", "GC4B", "GC4B", "GC5A", "GC5B", "GC5B",
# 	"LC5", "LC5L", "LC5L", "LC6", "LC6L", "LC6L", "LC11", "LC11L", "LC11L")
#
# control.list <- c("Breast", "Normal", "Blood", "GC1B", "GC1C", "GC1D", "GC2B", "GC2C", "GC2D", "GC3B", "GC3C", "GC3D", "GC4B", "GC4C", "GC4D", "GC5B", "GC5C", "GC5D",
# 	"LC5L", "LC5N", "LC5B", "LC6L", "LC6N", "LC6B", "LC11L", "LC11N", "LC11B")

# WGS
# sample.list <- c("D430P", "D430M", "D430M", "D441P", "D441M", "D441M", "D473P", "D473M", "D473M", "D559P", "D559M", "D559M")
# control.list <- c("D430N", "D430N", "D430P", "D441N", "D441N", "D441P", "D473N", "D473N", "D473P", "D559N", "D559N", "D559P")

# new breast samples
# sample.list <- c("Breast1_3", "Breast1_2", "Breast1_1", "Breast1_2", "Breast1_1", "Breast1_1", "Breast2_3", "Breast2_2", "Breast2_1", "Breast2_2", "Breast2_1", "Breast2_1", "Breast3_3", "Breast3_2", "Breast3_1", "Breast3_2", "Breast3_1", "Breast3_1")
# control.list <- c("Breast1_b", "Breast1_b", "Breast1_b", "Breast1_3", "Breast1_3", "Breast1_2", "Breast2_b", "Breast2_b", "Breast2_b", "Breast2_3", "Breast2_3", "Breast2_2", "Breast3_b", "Breast3_b", "Breast3_b", "Breast3_3", "Breast3_3", "Breast3_2")

# new lung samples
# sample.list <- c("LG39N", "LG39T", "LG39T")
# control.list <- c("LG39B", "LG39B", "LG39N")

# Breast06
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
