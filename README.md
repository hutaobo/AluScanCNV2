# AluScanCNV2: an R package for copy number variation-based cancer risk prediction

## Installation
AluScanCNV2 can be installed using the install_github function in the devtools package.
```{r}
library(devtools)
install_github('hutaobo/AluScanCNV')
```

## CNV calling
The coverageBed tool of the BEDtools software (Quinlan and Hall, 2010) was used to calculate the depth of sequencing reads in the analysis-ready BAM file (coverageBed -hist -a 5k.bin -b output.bed > output.5k.doc). The outputted coverage file was applied to the AluScanCNV software (Yang et al., 2014) for CNV calling.

AluScanCNV2 relies on Geary-Hinkley transformation (GHT)-based comparison of the read-depth of a sequence window on the test sample with that on either a paired control sample in the case of 'paired CNV' analysis, or a reference template constructed from pooled reference samples in the case of 'unpaired CNV' analysis (Yang et al., 2014).
```{r}
# Calling of paired CNV
library(AluScanCNV)
control_doc_path <- system.file("extdata/Breast1_b.5k.doc", package = "AluScanCNV")
tumor_doc_path <- system.file("extdata/Breast1_1.5k.doc", package = "AluScanCNV")
pairedCNV(control.5k.doc = control_doc_path, sample.5k.doc = tumor_doc_path, window.size = "500k", output.path = "./")

# Calling of unpaired CNV
sample_doc_path <- system.file("extdata/Breast1_b.5k.doc", package = "AluScanCNV")
unpairedCNV(sample.5k.doc = sample_doc_path, window.size = "500k", seq.method = "AluScan", output.path = "./")
```

## Identification of recurrent CNVs
```{r}
alu_control <- list.files(path = 'path_to_folder_of_seg_files', full.names = TRUE)
library(AluScanCNV)
alu_control <- seg2CNV(alu_control)
alu_control$recurrence <- alu_control$recurrence / (ncol(alu_control) - 4)
```

## Selection of informative CNVs
```{r}
recurr_cnv <- featureSelection2(nonCancerListA = alu_control, CancerListA = alu_cancer, nonCancerListB = wgs_control, CancerListB = wgs_cancer, Cri = 0.33)
```

## Prediction of cancer-predisposition
```{r}
library(caret)
metric <- "Accuracy"
control <- trainControl(method = "cv", number = 10)

# a) linear algorithms
fit.lda <- train(type ~ ., data = dataset, method = "lda", metric = metric, trControl = control)
# b) nonlinear algorithms
# CART
fit.cart <- train(type ~ ., data = dataset, method = "rpart", metric = metric, trControl = control)
# kNN
fit.knn <- train(type ~ ., data = dataset, method = "knn", metric = metric, trControl = control)
# c) advanced algorithms
# SVM
fit.svm <- train(type ~ ., data = dataset, method = "svmRadial", metric = metric, trControl = control)
# Random Forest
fit.rf <- train(type ~ ., data = dataset, method = "rf", metric = metric, trControl = control)
# d) others
# Naive Bayes
fit.nb <- train(type ~ ., data = dataset, method = "naive_bayes", metric = metric, trControl = control)

results <- resamples(list(lda = fit.lda, cart = fit.cart, knn = fit.knn, svm = fit.svm, rf = fit.rf, nb = fit.nb))
results <- resamples(list(cart = fit.cart, knn = fit.knn, svm = fit.svm, rf = fit.rf, nb = fit.nb))
summary(results)
```

## References
Quinlan, A. R. and I. M. Hall (2010). "BEDTools: a flexible suite of utilities for comparing genomic features." Bioinformatics 26(6): 841-842.

Yang, J. F., X. F. Ding, L. Chen, W. K. Mat, M. Z. Xu, J. F. Chen, J. M. Wang, L. Xu, W. S. Poon, A. Kwong, G. K. Leung, T. C. Tan, C. H. Yu, Y. B. Ke, X. Y. Xu, X. Y. Ke, R. C. Ma, J. C. Chan, W. Q. Wan, L. W. Zhang, Y. Kumar, S. Y. Tsang, S. Li, H. Y. Wang and H. Xue (2014). "Copy number variation analysis based on AluScan sequences." J. Clin. Bioinforma. 4(1): 15.
