# Methods-comparison-for-interomics-correlation-detection
The detection of microbial association is affected by community structure and sequencing depth. Some methods were specifically designed for the detection of associations within microbiome data. Currently, increasing attention is paid to the association between the microbiome and other omics data which requires strong supporting of appropriate methods. The R code and datasets used for performances comparison of 6 typical correlation methods were provided here, in view of the association detection on metabolome and microbiome data sets. 

Methods involved: Pearson, Spearman, SparCC (Sparse Correlations for Compositional data), CCLasso (Correlation inference for Compositional data through Lasso), MIC (Mutual Information Coefficient) and Cosine similarity.

Performances under comparison: specificity, sensitivity, similarity, accuracy, and stability on different sample sizes and zero value proportions. 

Some comments on methods selection could be concluded from the preliminary comparisons. 1) Each method has its own strength. Jointly usage of more than one method (especially those from different clusters) may achieve a better result. 2) Spearman is of optimal overall performance and maybe the first choice. MIC may serve as a good supporting method when the data set diversity, sample size, or the correlation strength fluctuates in a wide range. 3) Both zero values and sample size may affect the accuracy. The validity and reliability of result are low when there are more than 40% zero values or more than 50% of samples are missing. 4) Zero values have a stronger negative impact on correlation result compared with the reduction in sample size. Particular care should be taken when replacing missing values by zero values. 
 

