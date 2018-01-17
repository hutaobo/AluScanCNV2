#' Quantsmooth plot of selected CNV features
#' @param
#' @keywords
#' @export
#' @examples
#' plotSelectedFeatures()

plotSelectedFeatures <- function(output_path, nonCancerListA, CancerListA, nonCancerListB, CancerListB, Cri = 0.25) {
  library(quantsmooth)

  #- CancerListA <- list.files(path = '/mnt/data1/gene/thuac/Weka/result/seg_alu/tumor/500k/', full.names = TRUE)
  CancerListA <- seg2CNV(CancerListA)
  CancerListA$recurrence <- CancerListA$recurrence / (ncol(CancerListA) - 4)
  CancerListA_start <- data.frame(CHR = gsub("chr", "", CancerListA$chr), MapInfo = CancerListA$start, stringsAsFactors = FALSE)
  CancerListA_end <- data.frame(CHR = gsub("chr", "", CancerListA$chr), MapInfo = CancerListA$end, stringsAsFactors = FALSE)

  #- nonCancerListA <- list.files(path = '/mnt/data1/gene/thuac/Weka/result/seg_alu/control/500k/', full.names = TRUE)
  nonCancerListA <- seg2CNV(nonCancerListA)
  nonCancerListA$recurrence <- nonCancerListA$recurrence / (ncol(nonCancerListA) - 4)
  nonCancerListA_start <- data.frame(CHR = gsub("chr", "", nonCancerListA$chr), MapInfo = nonCancerListA$start, stringsAsFactors = FALSE)
  nonCancerListA_end <- data.frame(CHR = gsub("chr", "", nonCancerListA$chr), MapInfo = nonCancerListA$end, stringsAsFactors = FALSE)

  #- CancerListB <- list.files(path = '/mnt/data1/gene/thuac/Weka/result/seg/KG/500k/', full.names = TRUE)
  CancerListB <- seg2CNV(CancerListB)
  CancerListB$recurrence <- CancerListB$recurrence / (ncol(CancerListB) - 4)
  CancerListB_start <- data.frame(CHR = gsub("chr", "", CancerListB$chr), MapInfo = CancerListB$start, stringsAsFactors = FALSE)
  CancerListB_end <- data.frame(CHR = gsub("chr", "", CancerListB$chr), MapInfo = CancerListB$end, stringsAsFactors = FALSE)

  #- nonCancerListB <- list.files(path = '/mnt/data1/gene/thuac/Weka/result/seg/ICGC/500k/', full.names = TRUE)
  nonCancerListB <- seg2CNV(nonCancerListB)
  nonCancerListB$recurrence <- nonCancerListB$recurrence / (ncol(nonCancerListB) - 4)
  nonCancerListB_start <- data.frame(CHR = gsub("chr", "", nonCancerListB$chr), MapInfo = nonCancerListB$start, stringsAsFactors = FALSE)
  nonCancerListB_end <- data.frame(CHR = gsub("chr", "", nonCancerListB$chr), MapInfo = nonCancerListB$end, stringsAsFactors = FALSE)

  control <- merge(nonCancerListA, nonCancerListB, by = c("chr", "start", "end"), sort = FALSE)
  control$recurrence <- pmin(control$recurrence.x, control$recurrence.y)
  control_start <- data.frame(CHR = gsub("chr", "", control$chr), MapInfo = control$start, stringsAsFactors = FALSE)
  control_end <- data.frame(CHR = gsub("chr", "", control$chr), MapInfo = control$end, stringsAsFactors = FALSE)

  cancer <- merge(CancerListA, CancerListB, by = c("chr", "start", "end"), sort = FALSE)
  cancer$recurrence <- pmin(cancer$recurrence.x, cancer$recurrence.y)
  cancer_start <- data.frame(CHR = gsub("chr", "", cancer$chr), MapInfo = cancer$start, stringsAsFactors = FALSE)
  cancer_end <- data.frame(CHR = gsub("chr", "", cancer$chr), MapInfo = cancer$end, stringsAsFactors = FALSE)

  total <- merge(control, cancer, by = c("chr", "start", "end"), sort = FALSE)
  total <- total[, c('chr', 'start', 'end', 'recurrence.x.x', 'recurrence.y.x', 'recurrence.x', 'recurrence.x.y', 'recurrence.y.y', 'recurrence.y')]
  total$recurrence <- total$recurrence.x - total$recurrence.y
  total$color <- ifelse(abs(total$recurrence) > Cri, "red", "grey")
  total$up <- ifelse(total$recurrence > 0, abs(total$recurrence), 0)
  total$down <- ifelse(total$recurrence < 0, abs(total$recurrence), 0)
  total_start <- data.frame(CHR = gsub("chr", "", total$chr), MapInfo = total$start, stringsAsFactors = FALSE)
  total_end <- data.frame(CHR = gsub("chr", "", total$chr), MapInfo = total$end, stringsAsFactors = FALSE)

  threshold_start <- data.frame(CHR = 1:22, MapInfo = 1)
  threshold_end <- data.frame(CHR = 1:22, MapInfo = lengthChromosome(1:22))

  pdf(file = output_path, width = 7, height = 6)
  threshold_start_pos <- prepareGenomePlot(threshold_start, paintCytobands = FALSE, organism = "hsa", units = "hg19", sexChromosomes = FALSE)
  threshold_end_pos <- prepareGenomePlot(threshold_end, paintCytobands = FALSE, organism = "hsa", units = "hg19", sexChromosomes = FALSE)
  total_start_pos <- prepareGenomePlot(total_start, paintCytobands = FALSE, organism = "hsa", units = "hg19", sexChromosomes = FALSE)
  total_end_pos <- prepareGenomePlot(total_end, paintCytobands = TRUE, organism = "hsa", units = "hg19", sexChromosomes = FALSE)
  for (i in 1:nrow(total_start_pos)) {
    rect(total_start_pos[i, 2], total_start_pos[i, 1] + 0.06, total_end_pos[i, 2], total_end_pos[i, 1] + 0.06 + total$up[i] / 0.9, col = ifelse(total$color[i] == "grey", "steelblue", total$color[i]), border = NA)
    rect(total_start_pos[i, 2], total_start_pos[i, 1] - 0.06 - total$down[i] / 0.9, total_end_pos[i, 2], total_end_pos[i, 1] - 0.06, col = ifelse(total$color[i] == "grey", "steelblue", total$color[i]), border = NA)
  }
  for (i in 1:22) {
    lines(c(threshold_start_pos[i, 2], threshold_end_pos[i, 2]), c(threshold_start_pos[i, 1] + 0.06 + Cri / 0.9, threshold_end_pos[i, 1] + 0.06 + Cri / 0.9), lty = 2, col = "black", lwd = 0.5)
    lines(c(threshold_start_pos[i, 2], threshold_end_pos[i, 2]), c(threshold_start_pos[i, 1] - 0.06 - Cri / 0.9, threshold_end_pos[i, 1] - 0.06 - Cri / 0.9), lty = 2, col = "black", lwd = 0.5)
  }
  dev.off()
}
