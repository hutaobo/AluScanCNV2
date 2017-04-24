#' @param
#' @keywords internal
#' @export
#' @examples
#' z2t()

z2t <- function(z, lambdax, lambday) {
    (lambday * z - lambdax)/sqrt(lambday * z^2 + lambdax)
}


#' @param
#' @keywords internal
#' @export
#' @examples
#' cnv.ANNO()

cnv.ANNO <- function(cnv.sub, anno.last, step, minimum.window, threshold) {
    if (nrow(cnv.sub) > 1) {
        distance <- cnv.sub[2:nrow(cnv.sub), "start"] - cnv.sub[1:(nrow(cnv.sub) - 1), "start"]
        cnv.end <- c(which(distance > step), nrow(cnv.sub))
        last <- 1
        for (i in cnv.end) {
            if (i - last + 1 >= minimum.window) {
                anno.last <- anno.last + 1
                cnv.sub[last:i, "cnv"] <- anno.last
            }
            last <- i + 1
        }
    } else {
        if (nrow(cnv.sub) == 1) {
            if (minimum.window <= 1) {
                anno.last <- anno.last + 1
                cnv.sub$cnv <- anno.last
            }
        }
    }
    cnv.sub
}


#' https://github.com/hliang/cnv-seq/blob/master/cnv/R/cnv.R
#' @param
#' @keywords internal
#' @export
#' @examples
#' cnv.cal()

cnv.cal <- function(data, log2.threshold = 0.6, chromosomal.normalization = FALSE, annotate = FALSE, minimum.window = 4) {
    if (is.na(log2.threshold)) {
        log2.threshold <- as.numeric(sub("\\.pvalue-.+$", "", sub("^.+\\.log2-", "", file, perl = TRUE), perl = TRUE))
    }
    # data <- read.delim(file)
    data$position <- round((data$end + data$start)/2)
    data$log2 <- NA
    data$p.value <- NA
    test.lambda <- list()
    ref.lambda <- list()
    chrom <- as.character(unique(data$chromosome))
    if (chromosomal.normalization) {
        for (chr in chrom) {
            sub <- subset(data, chromosome == chr)
            lambda.test <- mean(sub$test)
            test.lambda[[chr]] <- lambda.test
            lambda.ref <- mean(sub$ref)
            ref.lambda[[chr]] <- lambda.ref
            norm <- lambda.test/lambda.ref
            ratio <- sub$test/sub$ref
            sub$log2 <- log2(sub$test/sub$ref/norm)
            data[rownames(sub), "log2"] <- sub$log2
            subp <- subset(sub, log2 >= 0)
            if (nrow(subp) > 0)
                data[rownames(subp), "p.value"] <- pnorm(z2t(subp$test/subp$ref, lambda.test, lambda.ref), lower.tail = FALSE)
            subn <- subset(sub, log2 < 0)
            if (nrow(subn) > 0)
                data[rownames(subn), "p.value"] <- pnorm(z2t(subn$test/subn$ref, lambda.test, lambda.ref), lower.tail = TRUE)
        }
    } else {
        # GClambdaTest has the same length with ratio
        lambda.test <- data$GClambdaTest
        lambda.ref <- data$GClambdaRef
        # lambda.test <- mean(data$test) lambda.ref <- mean(data$ref)
        for (chr in chrom) {
            # test.lambda[[chr]]<-lambda.test ref.lambda[[chr]]<-lambda.ref
            test.lambda[[chr]] <- mean(data$test)
            ref.lambda[[chr]] <- mean(data$ref)
        }
        # norm <- lambda.test/lambda.ref norm = mean(lambda.test)/mean(lambda.ref)
        norm <- sum(data$test)/sum(data$ref)
        ratio <- data$test/data$ref
        data$log2 <- log2(ratio/norm)
        subp <- subset(data, log2 >= 0)
        # data[rownames(subp),'p.value'] <- pnorm(z2t(subp$test/subp$ref,lambda.test,lambda.ref), lower.tail=FALSE)
        value <- z2t(subp$test/subp$ref, lambda.test[data$log2 >= 0], lambda.ref[data$log2 >= 0])
        data[rownames(subp), "p.value"] <- pnorm(value, sd = sd(value, na.rm = TRUE), lower.tail = FALSE)
        data[rownames(subp), "zScore"] <- value
        subn <- subset(data, log2 < 0)
        # data[rownames(subn),'p.value'] <- pnorm(z2t(subn$test/subn$ref,lambda.test,lambda.ref), lower.tail=TRUE)
        value <- z2t(subn$test/subn$ref, lambda.test[data$log2 < 0], lambda.ref[data$log2 < 0])
        data[rownames(subn), "p.value"] <- pnorm(value, sd = sd(value, na.rm = TRUE), lower.tail = TRUE)
        data[rownames(subn), "zScore"] <- value
    }
    if (annotate) {
        data$cnv <- 0
        data$cnv.size <- NA
        data$cnv.log2 <- NA
        data$cnv.p.value <- NA

        window.size <- data[1, "end"] - data[1, "start"] + 1
        step <- window.size/2
        for (chr in chrom) {
            print(paste("chromosome: ", chr))
            sub <- subset(data, chromosome == chr)
            cnv.p <- cnv.ANNO(subset(sub, log2 >= abs(log2.threshold)), max(sub$cnv, data$cnv), step, minimum.window, log2.threshold)
            if (nrow(subset(cnv.p, cnv > 0))) {
                data[rownames(cnv.p), "cnv"] <- cnv.p$cnv
            }
            cnv.n <- cnv.ANNO(subset(sub, log2 <= -1 * abs(log2.threshold)), max(sub$cnv, data$cnv), step, minimum.window, log2.threshold)
            if (nrow(subset(cnv.n, cnv > 0))) {
                data[rownames(cnv.n), "cnv"] <- cnv.n$cnv
            }
        }

        if (max(data$cnv) > 0) {
            for (id in seq(1, max(data$cnv))) {
                print(paste("cnv_id: ", id, " of ", max(data$cnv)))
                sub <- subset(data, cnv == id)
                start <- ceiling(mean(c(min(sub$start), min(sub$position))))
                end <- floor(mean(c(max(sub$end), max(sub$position))))
                size <- end - start + 1
                chr <- as.character(unique(sub$chromosome))
                ratio <- sum(sub$test)/sum(sub$ref)
                norm <- test.lambda[[chr]]/ref.lambda[[chr]]
                cnv.p.value <- 2 * pnorm(z2t(ratio, test.lambda[[chr]] * nrow(sub), ref.lambda[[chr]] * nrow(sub)), lower.tail = ratio < norm)
                index <- which(data$cnv == id)
                data[index, "cnv.log2"] <- log2(ratio/norm)
                data[index, "cnv.p.value"] <- cnv.p.value
                data[index, "cnv.size"] <- size
            }
        }
    }
    data
}
