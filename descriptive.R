# descriptive stat

# load data
load('stress_aligned.RData')

# raw response scores
# primary stressor
# secondary stressor
# Perceived Support
# group identity
# PSS
# BRS

library(psych)
data.filtered$STR1 <- rowMeans(data.filtered[,items.ps])
data.filtered$STR2 <- rowMeans(data.filtered[,items.ss])
data.filtered$SPS <- rowMeans(data.filtered[,items.sps])
data.filtered$PSS <- rowMeans(data.filtered[,items.pss])
data.filtered$BRS <- rowMeans(data.filtered[,items.resilience])
data.filtered$GID <- rowMeans(data.filtered[,items.identity])

vars <- c('STR1','STR2','SPS','PSS','BRS','GID')

bycountry <- describeBy(data.filtered[,vars],group=data.filtered$residing_country)

# create matrix
n.country <- length(table(data.filtered$residing_country))
result <- matrix(nrow=n.country, ncol = 13)
for (i in 1:n.country){
  for (j in 1:6){
    result[i,j*2] <- bycountry[[i]][j,3]
    result[i,j*2+1] <- bycountry[[i]][j,4]
  }
  
  
}
result[,1 ] <- as.numeric(table(data.filtered$residing_country))
result.table <- as.data.frame((result))
result.table <- cbind(labels(table(data.filtered$residing_country))[[1]],result.table)

# all countries
describe(data.filtered[,vars])

# alphas
psych::alpha(data.filtered[,items.ps],check.keys = T)
psych::alpha(data.filtered[,items.ss],check.keys = T)
psych::alpha(data.filtered[,items.sps],check.keys = T)
psych::alpha(data.filtered[,items.identity],check.keys = T)
psych::alpha(data.filtered[,items.pss],check.keys = T)
psych::alpha(data.filtered[,items.resilience],check.keys = T)
