# Returns the average fragment length flanked by MseI and EcoRI, for genome size G and GC-content GC.
avgFragLength <- function(G, GC){
  library(Digestion)
  library(dplyr)
  data(allRestrictionEnzymes)
  genome <- paste(sample(x=c("A","C","G","T"),                 
                         size=G,                               
                         replace=T,                            
                         prob=c((1-GC)/2,GC/2,GC/2,(1-GC)/2)),
                  collapse="")
  fragments <- getDigestionFragments(Restriction=c("EcoRI:MseI"), DNASequence=genome)
  fragments_1 <- filter(fragments, (fragments$EnzymeStart == "MseI" & fragments$EnzymeEnd == "EcoRI") | (fragments$EnzymeStart == "EcoRI" & fragments$EnzymeEnd == "MseI"))
  return(mean(fragments_1$Length))
}


test <- vector(mode="numeric", length=0)
for(i in 1:10000) {
  test[i] <- avgFragLength(1000000, 0.5)
}
mean(test) # Estimated fragment length


h<-hist(test, breaks=200, col="red", xlab="Fragment length", 
        main="Histogram with Normal Curve") 
xfit<-seq(min(test),max(test),length=40) 
yfit<-dnorm(xfit,mean=mean(test),sd=sd(test)) 
yfit <- yfit*diff(h$mids[1:2])*length(test) 
lines(xfit, yfit, col="blue", lwd=2)