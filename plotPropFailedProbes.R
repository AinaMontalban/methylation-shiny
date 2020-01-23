
plotFailedPropProbes<- function(detP, sampleNames){
  pal <- brewer.pal(8,"Paired")
  
failed_proportion <- apply(detP, 2, function(x) {
  # All probes detP == 0 means totally failed array.
  if (sum(x) == 0) {
    return(NA)
  } else {
    failed_probes <- sum(x >= 0.01)
    return(failed_probes / length(x))
  }
})

names(failed_proportion) <- sampleNames
barplot(failed_proportion, col=pal[factor(sampleNames)], 
        ylim = c(0, 0.12),
        las=2, 
        cex.names=0.8, ylab="Proportion of failed probes")
abline(h=0.1,col="red")
abline(h=0.05,col="blue")
abline(h= 0.01, col = "green")
}
