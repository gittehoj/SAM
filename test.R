library(stockassessment)
example(sam.fit)

conf<-nscodConf
conf$stockWeightProcess <- 2
par <- defpar(nscodData, conf)
fit2 <- sam.fit(nscodData, conf, par)

logSW<-fit2$pl$logSW
for(i in 1:ncol(logSW)){
  for(j in 1:nrow(logSW)){
    if(i>=j)logSW[j,i]<-logSW[j,i]+fit2$pl$logCSW[i-j+1]
  }
}

matplot(nscodData$stockMeanWeight)
matplot(t(exp(logSW)), type="l", add=TRUE)
