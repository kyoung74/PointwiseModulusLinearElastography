mapFitToTXT=function(mapFits){
  # Takes the object produced by parExtractStiffness in Rasylum
  noFits=length(mapFits$fits)
  for(f in 1:noFits){
    print(paste(f," of ",noFits,sep=""))
    truncated=mapFits$fits[[f]]$fit$curves[mapFits$fits[[f]]$fit$curves$curve=="measured"&mapFits$fits[[f]]$fit$curves$zPos>=0,]
    write.table(cbind(truncated$zPos,truncated$F),file=paste("Line",mapFits$fits[[f]]$ident$Line,"Point",mapFits$fits[[f]]$ident$Point,".txt",sep=""),sep="\t",row.names=FALSE,col.names=c("Indentation","Force"))
  }
}