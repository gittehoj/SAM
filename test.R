library(stockassessment)
conf<-nscodConf
conf$stockWeightProcess <- 1
par <- defpar(nscodData, conf)
fit <- sam.fit(nscodData, conf, par)
matplot(nscodData$stockMeanWeight)
matplot(t(exp(fit$pl$logSW)), type="l", add=TRUE)
