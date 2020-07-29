## Convert beta-values to M-values
beta <- readRDS("../data/raw_data/beta_cleaned_NL.rds")

onevalues <- which(beta == 1)
zerovalues<- which(beta == 0)

if(length(onevalues) > 0 || length(zerovalues) >0 ) {
  
  if(length(onevalues)>0) {
    beta[onevalues] <- NA
    beta[onevalues] <- max(beta, na.rm=T)
  }
  if(length(zerovalues) > 0)  {
    beta[zerovalues] <- NA
    beta[zerovalues] <- min(beta, na.rm=T)
  }
}

mvalues = log2(beta/(1-beta))

saveRDS(mvalues, file = "../data/raw_data/mvalues_cleaned_NL.rds")