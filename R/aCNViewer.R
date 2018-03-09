#' Plot the frequency of CNV
#' @param
#' @keywords
#' @export
#' @examples
#' plotFrequency()

plotFrequency <- function(input, Cri = 0.05, return = TRUE, output_path = NULL, recCri = NULL) {
  require(GenomicRanges)
  require(data.table)
  require(ggplot2)
  cnv_gr <- c()
  for (i in 1:length(input)) {
    data <- fread(input[i], header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
    cnv <- data.frame(chr = data$chromosome, start = data$start, end = data$end, cnv = ifelse(data$p.value >= Cri, 0, ifelse(data$zScore > 0, 1, -1)), stringsAsFactors = FALSE)
    cnv <- subset(cnv, cnv != 0)
    cnv_gr <- rbind(cnv_gr, cnv)
  }
  cnv_gr <- makeGRangesFromDataFrame(cnv_gr, keep.extra.columns = TRUE)

  loss_gr <- subset(cnv_gr, cnv == -1)
  loss_cov <- disjoin(loss_gr)
  loss_cov$coverage <- countOverlaps(loss_cov, loss_gr) * -1

  gain_gr <- subset(cnv_gr, cnv == 1)
  gain_cov <- disjoin(gain_gr)
  gain_cov$coverage <- countOverlaps(gain_cov, gain_gr)

  cnv <- c(gain_cov, loss_cov)

  pdata <-
    data.frame(
      chr = gsub("chr", "", c(seqnames(cnv))),
      xmin = start(cnv),
      xmax = end(cnv),
      height = cnv$coverage / length(input),
      fill = ifelse(cnv$coverage > 0, "gain", "loss")
    )

  if (is.null(recCri)) {
    recCri <- 0.1
  }
  p <- ggplot(pdata, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = height, fill = fill)) +
    geom_rect() +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
    theme(panel.spacing = unit(0.1, "lines")) +
    geom_hline(aes(yintercept = recCri), linetype = "dashed") +
    geom_hline(aes(yintercept = -recCri), linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 30)) +
    scale_color_manual(c("gain" = "green", "loss" = "orange"))

  if (return) {
    return(p)
  } else {
    pdf(file = file.path(output_path, "Frequency_Of_CNV.pdf"))
    plot(p)
    dev.off()
  }
}
