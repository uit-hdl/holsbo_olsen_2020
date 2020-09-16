library(rstan)

options(mc.cores = 10)

load("data/metastasis.rda")

# transform all genes to be on same scale
fchange <- t(plyr::aaply(fchange, 2, function(x) (x-mean(x))/sd(x)))

# change from factor type to integers (for stan)
covariables$status <- as.factor(as.character(covariables$status))
metastasis <- ifelse(covariables$metastasis=="Yes", 1, 0)
detection <- ifelse(covariables$status == "screening", 1, 2)

expression <- t(fchange)
s_data <- list(expression=expression, metastasis=metastasis, N=ncol(expression), 
               M=nrow(expression), L=2, detection=detection)

# a check
mns <- -(rowMeans(expression[, metastasis==1]) - rowMeans(expression[, metastasis==0]))
hist(mns)
quantile(mns)
##           0%          25%          50%          75%         100% 
## -0.608677155 -0.092230287  0.009288161  0.103710246  0.504351198 


#        90% 
# 0.03545651 

# subset of data used to see that the model runs for fewer genes, not input into final model
s_data_test <- s_data
s_data_test$expression <- s_data$expression[1:1000, ]
s_data_test$M <- 1000

# compile and run model, save results to file
stanc("models/genewise_logistic.stan")
fit <- stan(
  "models/genewise_logistic.stan",
  data = s_data,
  chains = 5,
  warmup = 1000,
  iter = 2200,
  
); save(fit, file = "data/genewise_logistic.rda")

 

