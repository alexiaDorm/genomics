library("survival")
library(survminer)
library(ggfortify)

univariate_Cox <- function(test_data, genes){
  ranked_genes <- vector(mode="double", length=length(genes)); names(ranked_genes) <- genes
  for (i in 1:length(genes)){
    formula <- as.formula(paste('Surv(DMFSVal, DMFS, type="right")~', genes[i]))
    res.cox <- coxph(formula, data = test_data)
    x <- summary(res.cox)
    ranked_genes[i] <- x$wald["pvalue"]
  }
  ranked_genes <- sort(ranked_genes, decreasing = F)
  
  return (ranked_genes)
}

compute_DMFS <- function(exp){
  #Compute percentage DMFS for each group
  DMFS.low <- (sum(exp[exp$group == 1,]$DMFS)/length(exp[exp$group == 1,]$DMFS))
  DMFS.medium <- (sum(exp[exp$group == 2,]$DMFS)/length(exp[exp$group == 2,]$DMFS))
  DMFS.high <- (sum(exp[exp$group == 3,]$DMFS)/length(exp[exp$group == 3,]$DMFS))
  
  DMFS.groups <- c(DMFS.low,DMFS.medium,DMFS.high)
  DMFS.groups <- 1- DMFS.groups 
  DMFS.groups[is.na(DMFS.groups)] <- 0
  
  return (DMFS.groups)
}

assign_exp_group<- function(exp){
  #Assign low,medium,high to each sample
  exp$sum <- rowMeans(exp)
  tert = quantile(exp$sum, c(0:3/3))
  exp$group = with(exp, cut(sum, tert, include.lowest = T, labels = c("Low", "Medium", "High")))
  
  return (exp$group)
}

refine_metagene<- function(exp, genes){
  #Shuffle the data and create fold
  exp <- exp[sample(nrow(exp)),]
  folds <- cut(seq(1,nrow(exp)),breaks=10,labels=FALSE)
  
  #Perform 10 fold CV
  exp_group <- data.frame(matrix(0, nrow = nrow(exp), ncol = length(genes)-9))
  colnames(exp_group) <- 10:length(genes); rownames(exp_group) <- rownames(exp)
  for(i in 1:10){
    #Get train and test data sets
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- exp[testIndexes, ]
    trainData <- exp[-testIndexes, ]
    
    #Perform univariate Cox analysis of each gene on train set
    ranked_genes <- univariate_Cox(trainData, genes)
    
    #Assign group to test sets for each metagene size
    for(N in 10:length(genes)){
      metagene <- names(ranked_genes)[1:N]
      exp_group[testIndexes, N-9] <- assign_exp_group(testData[, metagene])
    }
  }
  
  #Compute 5 years DMFS (%) and log-rank test for each size of metagene
  DMFS_metagenes = c(); p_value = c()
  for(N in 10:length(genes)){
    data <- data.frame(exp_group[N-9], exp$DMFS, exp$DMFSVal, row.names = rownames(exp))
    colnames(data) <- c("group","DMFS", "time")
    
    #Compute 5 years DMFS (%)
    DMFS_metagenes <- append(DMFS_metagenes, compute_DMFS(data))
    
    #Log-rank test
    surv_diff <- survdiff(Surv(time, DMFS)~ group, data = data)
    p_value <- append(p_value, (1 - pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = TRUE)))
  }
  
  return(c(DMFS_metagenes, p_value))
}