#################### FUNCTIONS for plotting megamate allele frequency data ################
#last edited 20160630 splitting G from smG functions


cumPosition<-function(myData,genomeIdx){
  cumPos=with(myData, data.frame("Chr"=Chr,"Pos"=Pos,"Position"=rep(0,dim(myData)[1])))
  levels(cumPos$Chr)<-c(levels(cumPos$Chr),levels(genomeIdx[,1])[5])
  genomeIdx<-cbind(genomeIdx,V6=cumsum(c(0,genomeIdx$V2))[1:dim(genomeIdx)[1]])
  for (i in 1:dim(genomeIdx)[1]){
    chr<-genomeIdx[i,1]
    cumPos[cumPos$Chr==chr,3]<-cumPos$Pos[cumPos$Chr==chr]+genomeIdx[genomeIdx$V1==chr,6]
  }
  myData<-cbind(myData[,1:3],cumPos[3],myData[,4:dim(myData)[2]])
  return(myData)
}

######### functions to find positions of genes of interest for plotting
cumPosition_mini<-function(Chr,Pos,genomeIdx){
  cumPos=data.frame("Chr"=Chr,"Pos"=Pos,"Position"=rep(0,length(Pos)))
  levels(cumPos$Chr)<-c(levels(cumPos$Chr),setdiff(levels(genomeIdx[,1]),levels(cumPos$Chr)))
  genomeIdx<-cbind(genomeIdx,V6=cumsum(c(0,genomeIdx$V2))[1:dim(genomeIdx)[1]])
  for (i in 1:dim(genomeIdx)[1]){
    chr<-genomeIdx[i,1]
    cumPos[cumPos$Chr==chr,3]<-cumPos$Pos[cumPos$Chr==chr]+genomeIdx[genomeIdx$V1==chr,6]
  }
  return(cumPos)
}

getLocNames<-function(gr,phenotypes,phenotype=NULL) {
  return(as.character(gr[gr$public %in% phenotypes[[phenotype]],]$public))
}

getLocPositions<-function(gr,phenotypes,phenotype,genomeIdx) {
  pos<-as.numeric(gr[gr$public %in% phenotypes[[phenotype]],]$locPos)
  chr<-as.character(seqnames(gr[gr$public %in% phenotypes[[phenotype]],]))
  return(cumPosition_mini(chr,pos,genomeIdx)$Position)
}
##################

cleanData<-function(myData,dataName,genomeVer,MADs=5,minReads=2,DelNonPolymorphic=FALSE) {
  myMessage<-paste(dataName,Sys.time(),"\n")
  #write to log file some general data
  myMessage<-paste0(myMessage,"All loci: ",dim(myData)[1],"\n")
  readDepth<-grep("_readDepth",names(myData))
  med<-median(myData[,readDepth])
  myMessage<-paste0(myMessage,"Median read depth: ",med," \n")
  mad<-mad(myData[,readDepth])
  myMessage<-paste0(myMessage,"Median average deviation of read depth: ",mad," \n")
  #remove sequences that are extreme outliers in the distribution
  remove<-myData[,readDepth]>med+MADs*mad
  myData<-myData[!remove,]
  myMessage<-paste0(myMessage,sum(remove)," loci removed because more than ", med+MADs*mad," reads (",MADs,"xMADs)\n")
  remove<-myData[,readDepth]<med-MADs*mad
  myData<-myData[!remove,]
  myMessage<-paste0(myMessage,sum(remove)," loci removed because less than ", med-MADs*mad," reads (",MADs,"xMADs)\n")
  #remove sequences with readDepth < minReads
  remove<-myData[,readDepth]<minReads
  myData<-myData[!remove,]
  myMessage<-paste0(myMessage,sum(remove)," loci removed because readDepth<",minReads,"\n")
  #remove loci that are not polymorphic
  if (DelNonPolymorphic) {
    CBfreq<-grep("_CBfreq",names(myData))
    remove<-(myData[,CBfreq] == 1 | myData[,CBfreq] == 0)
    myData<-myData[!remove,]
    myMessage<-paste0(myMessage,sum(remove)," loci removed because non-polymorphic in ", dataName,"\n")
  }
  #write some final data to logFile
  myMessage<-paste0(myMessage,"Final number of loci: ",dim(myData)[1],"\n")
  logFile<-file(description=paste0("../finalData/",genomeVer,"/log_",dataName),open="a")
  writeLines(myMessage,logFile)
  close(logFile)
  return(myData)
}



removeNonPolymorphic<-function(myData,dataName,genomeVer) {
  #remove sequences that are not polymorphic
  CBfreq<-grep("_CBfreq",names(myData))
  remove<-(rowMeans(myData[,CBfreq]==1)==1 | rowMeans(myData[,CBfreq]==0)==1)
  myData<-myData[!remove,]
  #write info to logfile
  myMessage<-paste(dataName,Sys.time(),"\n")
  myMessage<-paste0(myMessage,sum(remove)," loci removed because non-polymorphic in ", dataName,"\n")
  myMessage<-paste0(myMessage,"Final number of loci: ",dim(myData)[1],"\n")
  logFile<-file(description=paste0("../finalData/",genomeVer,"/log_",dataName),open="a")
  writeLines(myMessage,logFile)
  close(logFile)
  return(myData)
}



plotFreq <- function (myData, extraData, useLoci=c(2)) {
  CBfreq<-grep("_CBfreq",names(myData))
  plot(myData$Position,myData[,CBfreq[1]], type='p',xlab="Position(bp)",ylab="Hawaii allele frequency",
       main=extraData$mainTitle,ylim=c(0,1),xlim=c(0,max(myData$Position)),pch=16,cex=0.4,col="#00000077")
  points(myData$Position,myData[,CBfreq[2]],type='p',pch=16,cex=0.4,col="#FF000077")
  text(extraData$midpoint,0.98,extraData$chrNum,cex=1, col="black")
  abline(v=extraData$chrEnds,col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="dark red",lwd=0.6)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=0.7, col="red",
        adj=extraData$labelAdj[useLoci],las=2)
  legend("bottomright",legend=c("1F3","mel-26"),col=c(1,2),pch=16,cex=0.9)
}




plotSmoothedDiff <- function (myData, extraData, sm=1001, useLoci=c(2)) {
  CBfreq<-grep("_CBfreq",names(myData))
  CBdiff<-myData[,CBfreq[2]]-myData[,CBfreq[1]]
  smdiff<-smootheByChr(CBdiff,myData$Chr,sm)
  plot(myData$Position,smdiff, type="n",xlab="Position(bp)", ylab="Hawaii allele freq difference", 
       main=extraData$mainTitle, ylim=c(-1,1), xlim=c(0,max(myData$Position)))
  plotByChr(myData, smdiff, chrList=extraData$chrs,lwd=4)
  text(extraData$midpoint,0.98,extraData$chrNum,cex=1, col="black")
  abline(v=extraData$chrEnds,col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=0.7, col="red",
        adj=extraData$labelAdj[useLoci],las=2)
  abline(h=0,lty=5,col="dark gray")
}


smootheByChr <- function( data2smoothe, chrData, smWin) {
  smoothedData<-c(sgolayfilt(data2smoothe[chrData=="CHROMOSOME_I"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_II"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_III"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_IV"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_V"],p=3,n=smWin),
                  sgolayfilt(data2smoothe[chrData=="CHROMOSOME_X"],p=3,n=smWin))
  return(smoothedData)
}


############## old functions ######################
# values are not independant in adjacent position, so this is NOT a correct smoothing 
# Method
# apply fisher method with a rolling window. Adapted from MADAM package function. requires zoo package.
# adds "padding" values (NAs) to maintain original length of vector

## **** need to correct this. Since i am doing it one chromosome at a time, my k is too small!!!!!!!!
rollingFisher <- function (pvals, window=9, zero.sub = 1e-30, na.rm = FALSE) {
    stopifnot(all(pvals >= 0, na.rm = TRUE) & all(pvals <= 1, na.rm = TRUE))
    stopifnot(zero.sub >= 0 & zero.sub <= 1 || length(zero.sub) != 1)
    suppressPackageStartupMessages(library(zoo))
    pvals[pvals == 0] <- zero.sub
    logPs<- log(pvals)
    S=rollapply(logPs,window,function(x) {-2*sum(x)},partial=TRUE, align="center")
    fisher.sums<-data.frame(S=S,num.p=window)
    fisher.sums$p.value <- 1 - pchisq(fisher.sums$S, df = 2 *
                                          fisher.sums$num.p)
    return(fisher.sums)
}

# combine p values along chromosome using fisher method
fisherMethodByChr<-function(pValues, myData, sm=9, p.corr="boneferroni") {
    chrs<-levels(myData$Chr)
    FisherResults<-rollingFisher(pValues[myData$Chr==chrs[1]], sm)
    for (j in chrs[2:6]) {
        FisherResults<-rbind(FisherResults,rollingFisher(pValues[myData$Chr==j], sm))
    }
    return(p.adjust(FisherResults$p.value,method=p.corr))
}

plotPvals<-function(myData, log10pVals, extraData, useLoci=c(2)) {
    plot(myData$Position ,log10pVals, type="n",xlab="Position(bp)", ylab="-log10(p value)", main=extraData$mainTitle, 
         ylim=c(0,max(log10pVals)), xlim=c(0,max(myData$Position)))
    plotByChr(myData, log10pVals, chrList=extraData$chrs, lwd=4)
    text(extraData$midpoint, 0.98*max(log10pVals), extraData$chrNum, cex=1, col="black")
    abline(v=genomeIdx[,3],col="dark gray",lty=5)
    abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
    mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=1, col="red",adj=extraData$labelAdj[useLoci])
    abline(h=-log10(0.05),lty=5,col="dark red")
    text(90000000,-log10(0.05),"FDR=0.05", col="dark red",pos=3)
}
############### /old functions ###################



findPeaksByChr_slope <- function( data2smoothe, chrName, Position, smWin, slWin, weighted) {
  chr1<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_I"],Position[chrName=="CHROMOSOME_I"],
                        smWin,slWin,weighted)
  chr2<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_II"],Position[chrName=="CHROMOSOME_II"],
                        smWin,slWin,weighted)
  chr3<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_III"],Position[chrName=="CHROMOSOME_III"],
                        smWin,slWin,weighted)
  chr4<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_IV"],Position[chrName=="CHROMOSOME_IV"],
                        smWin,slWin,weighted)
  chr5<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_V"],Position[chrName=="CHROMOSOME_V"],
                        smWin,slWin,weighted)
  chrX<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_X"],Position[chrName=="CHROMOSOME_X"],
                        smWin,slWin,weighted)
  peakVals<-cbind(chr1[[1]],chr2[[1]],chr3[[1]],chr4[[1]],chr5[[1]],chrX[[1]])
  return(peakVals)
}


findPeaks_slope<-function(myVals,Position,smWin,slopeWin,weighted=TRUE) {
  #first smoothe the values with a moving window (smWin)
  smVals<-rollapply(myVals,smWin,mean)
  smVals<-c(rep(0,smWin/2),smVals,rep(0,smWin/2-1))
  #add weights (optional)
  if (weighted) {
    weights<-c(1:(slopeWin/2),rev(1:(slopeWin/2)))/(slopeWin/2)
  } else {
    weights<-NULL
  }
  #calculate slopes with a moving window (slopeWin)
  slopes<-rep(0,length(smVals)-slopeWin)
  for (i in 1:length(slopes)) {
    slopes[i]<-lm(smVals[i:(i+slopeWin-1)]~Position[i:(i+slopeWin-1)], weights=weights)$coefficients[2]
  }
  slopes<-c(rep(0,slopeWin/2),slopes,rep(0,slopeWin/2))
  #find slope sign
  slopeSign<-slopes>=0
  #find where the slope sign changes from + to - (i.e. diff: FALSE-TRUE= 0-1 = -1)
  peaks<-which(c(0,diff(slopeSign))==-1)
  return(list(rbind(smVals[peaks],Position[peaks]),smVals))
}


findPeaksByChr_slope_plot <- function( data2smoothe, chrName, Position, smWin, slWin, weighted,extraData=NULL) {
  chr1<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_I"],Position[chrName=="CHROMOSOME_I"],
                        smWin,slWin,weighted)
  chr2<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_II"],Position[chrName=="CHROMOSOME_II"],
                        smWin,slWin,weighted)
  chr3<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_III"],Position[chrName=="CHROMOSOME_III"],
                        smWin,slWin,weighted)
  chr4<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_IV"],Position[chrName=="CHROMOSOME_IV"],
                        smWin,slWin,weighted)
  chr5<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_V"],Position[chrName=="CHROMOSOME_V"],
                        smWin,slWin,weighted)
  chrX<-findPeaks_slope(data2smoothe[chrName=="CHROMOSOME_X"],Position[chrName=="CHROMOSOME_X"],
                        smWin,slWin,weighted)
  peakVals<-cbind(chr1[[1]],chr2[[1]],chr3[[1]],chr4[[1]],chr5[[1]],chrX[[1]])
  smVals<-c(chr1[[2]],chr2[[2]],chr3[[2]],chr4[[2]],chr5[[2]],chrX[[2]])
  #return(peakVals,smVals)
  myData<-data.frame(Chr=chrName,Position=Position)
  chrList<-levels(myExp$Chr)
  plot(Position,smVals,type='n',lwd=2,ylab="Smoothed -log10(p value)",main=extraData$mainTitle)
  plotByChr(myData, smVals, chrList, chrColumn="Chr",lwd=2)
  #plot(Position,smVals,type='l',lwd=2)
  #title(paste0("Smoothing Win=",smWin," slope Win=",slWin," weighted=",weighted))
  #abline(v=c(4185062,4189761),col="red")
  #points(peakVals[2,],peakVals[1,],col="black",pch=16,cex=0.5)
  text(extraData$midpoint,0,extraData$chrNum,cex=1, col="black")
  abline(v=extraData$chrEnds,col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=0.7, col="red",
        adj=extraData$labelAdj[useLoci])
  print(peakVals)
  return(min(abs(peakVals[2,]-mean(c(4185062,4189761)))))
}


plotByChr<-function(myData, yData, chrList, chrColumn="Chr", ...) {
  colours<-brewer.pal(6,"Dark2")
  for (c in chrList) {
    i=which(myData[,chrColumn]==c)
    lines(myData[i,"Position"], yData[i], col=colours[match(c,chrList)],...)
  }
}



plotSmPvals<-function(myData, log10pVals, extraData, sm=1001, useLoci=c(2)) {
  smPvals<-smootheByChr(log10pVals, myData$Chr, sm)
  plot(myData$Position ,smPvals, type="n",xlab="Position(bp)", ylab="-log10(p value)", 
       main=extraData$mainTitle, ylim=c(0,max(smPvals)), xlim=c(0,max(myData$Position)))
  plotByChr(myData, smPvals, chrList=extraData$chrs, lwd=4)
  text(extraData$midpoint, 0.98*max(smPvals), extraData$chrNum, cex=1, col="black")
  abline(v=extraData$chrEnds,col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=0.7, col="red",
        adj=extraData$labelAdj[useLoci],las=2)
  abline(h=-log10(0.05),lty=5,col="dark red")
  text(90000000,-log10(0.05),"FDR=0.05", col="dark red",pos=3)
}

plotSmGvals_sGolay<-function(myData, Gvals, extraData, sm=1001, useLoci=c(2)) {
  smGvals<-smootheByChr(Gvals, myData$Chr, sm)
  plot(myData$Position ,smGvals, type="n",xlab="Position(bp)", ylab="smoothed G values (sGolay)", 
       main=extraData$mainTitle, ylim=c(0,max(smGvals)), xlim=c(0,max(myData$Position)))
  plotByChr(myData, smGvals, chrList=extraData$chrs, lwd=4)
  text(extraData$midpoint, 0.98*max(smGvals), extraData$chrNum, cex=1, col="black")
  abline(v=extraData$chrEnds,col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=0.7, col="red",
        adj=extraData$labelAdj[useLoci],las=2)
}

plotVals<-function(myData, yvals, distanceMeasure, extraData, useLoci=c(2), fdr=NULL) {
  chrList<-levels(myData$Chr)
  plot(myData$Position,yvals,type='n',lwd=2,ylab="Smoothed G value",
     xlab="Position (bp)",main=extraData$mainTitle,ylim=c(0,max(yvals)*1.1))
  plotByChr(myData, yvals, chrList, chrColumn="Chr",lwd=2)
  text(extraData$midpoint,max(yvals)*1.1,extraData$chrNum,cex=1, col="black")
  abline(v=extraData$chrEnds,col="dark gray",lty=5)
  abline(v=extraData$locs[useLoci], col="red",lwd=0.8)
  mtext(extraData$locNames[useLoci], at=extraData$locs[useLoci], cex=0.7, col="red",
      adj=extraData$labelAdj[useLoci])
  title(sub=distanceMeasure)
}



library(DescTools, quietly=TRUE)
doGTest<-function(vector,...) {
  return(GTest(matrix(unlist(vector),nrow=2),...)$statistic)
}

getGbyChr<-function(myData,position,...) {
  library(DescTools,quietly=TRUE)
  chrs<-levels(myData$Chr)
  countCols<-grep("counts",names(myData))
  allGvals<-c()
  for (c in chrs) {
    gvals<-apply(myData[myData$Chr==c,countCols],1,doGTest,correct="williams")
    allGvals<-c(allGvals,gvals)
  }
  return(allGvals)
}

smootheGbyChr<-function(myData,winSize,position,gvalCol,...) {
  chrs<-levels(myData$Chr)
  allSmG<-c()
  for (c in chrs) {
    Gs<-myData[myData$Chr==c,gvalCol]
    pos<-position[myData$Chr==c]
    if (length(winSize)>1){
      win<-winSize[myData$Chr==c]
    } else {
      win<-winSize
    }
    smG<-smootheGvals(Gs,pos,win,...)
    allSmG<-c(allSmG,smG)
  }
  return(allSmG)
}

getG_smootheGbyChr<-function(myData,winSize,position,...) {
  library(DescTools,quietly=TRUE)
  chrs<-levels(myData$Chr)
  countCols<-grep("counts",names(myData))
  allSmG<-c()
  allGvals<-c()
  for (c in chrs) {
    gvals<-apply(myData[myData$Chr==c,countCols],1,doGTest,correct="williams")
    pos<-position[myData$Chr==c]
    if (length(winSize)>1){
      win<-winSize[myData$Chr==c]
    } else {
      win<-winSize
    }
    smG<-smootheGvals(gvals,pos,win,...)
    allGvals<-c(allGvals,gvals)
    allSmG<-c(allSmG,smG)
  }
  myData<-cbind(myData,Gval=allGvals,Gprime=allSmG)
  return(myData)
}


rr101byChr<-function(myData,recRate,halfWinSize,...) {
  library(DescTools,quietly=TRUE)
  chrs<-levels(myData$Chr)
  #countCols<-grep("counts",names(myData))
  rr101<-c()
  #allGvals<-c()
  for (c in chrs) {
    position<-1:(length(myData$Position[myData$Chr==c])+2*halfWinSize)
    rr<-c(rep(0,halfWinSize),recRate[myData$Chr==c],rep(0,halfWinSize))
    for (i in (halfWinSize+1):(length(position)-halfWinSize)) {
      rr101<-c(rr101,mean(rr[(i-halfWinSize):(i+halfWinSize)]))
    }
  }
  return(rr101)
}


tricubekernel<-function(distances) {
  Sw=sum((1-distances^3)^3)
  k=((1-distances^3)^3)/Sw
  return(k)
}

triweightkernel<-function(distances) {
  Sw=sum((1-distances^2)^3)
  k=((1-distances^2)^3)/Sw
  return(k)
}

biweightkernel<-function(distances) {
  Sw=sum((1-distances^2)^2)
  k=((1-distances^2)^2)/Sw
  return(k)
}

gaussiankernel<-function(distances) {
  Sw=sum((1/sqrt(2*pi)*exp(-0.5*distances^2)))
  k=(1/sqrt(2*pi)*exp(-0.5*distances^2))/Sw
  return(k)
}

smootheGvals<-function(gVals,pos,winSize,kernel=tricubekernel) {
  Gprime=rep(0,length(gVals))
  for (i in 1:length(gVals)){
    if (length(winSize)>1){
      j<-i
    } else {
      j=1
    }
    D=abs(pos-pos[i]) #calculate distance of all positions from current position
    winIndex=D<winSize[j]/2 #get index for all positions within 1/2 window size distance from current position
    Dwin=D[winIndex]*2/winSize[j] # standardize distances by window size (so values vary beween 0 and 1)
    k=kernel(Dwin) #get kernel weighting for each position in window
    Gprime[i]=sum(k*gVals[winIndex]) #calculate the weighted Gprime value by multiplying Gvals by kernel weightings
  }
  return(Gprime)
}

findDistppw1<-function(myData,ppw1_loc=c(4185062,4189761)) {
  ind<-which.max(myData$Gprime[myData$Chr=="CHROMOSOME_I"])
  dist1<-round(min(abs(myData$Position[ind]-ppw1_loc))/1000,1)
  put<-data.frame(x=ppw1_loc[2],y=max(myData$Gprime[myData$Chr=="CHROMOSOME_I"]))
  text(put$x,put$y,paste("max of peak is",dist1, "kb from ppw-1"), col="blue",cex=0.8,pos=4)
}


predictRecDist<-function(myData,RRsplineObj) {
  # Need the following two lines in script before calling this function:
  # load("RecRateSplitObj.RData") #loads object called s which is the smoothed spline of recombination rate
  # names(s)<-chrs #rename chromosomes in list
  recFrac<-c()
  chrs<-names(RRsplineObj)
  for (i in 1:length(RRsplineObj)) {
    b<-predict(RRsplineObj[[i]],myData$Pos[myData$Chr==chrs[i]])
    delta<-c(0,diff(b$y)) #find all positions where b$y is non-increasing 
    delta[delta<0]<-0 # replace negative increments with 0
    bnew<-cumsum(delta) # regenerate cumulative recombination fraction
    recFrac<-c(recFrac,bnew)
  }
  return(recFrac)
}

binData<-function(myData,binFactor,binSize) {
  countCols<-grep("counts",names(myData))
  newTable<-c()
  for (c in levels(myData$Chr)) {
    thisChr<-myData[myData$Chr==c,]
    bF<-binFactor[myData$Chr==c]
    binBorders<-pretty(c(min(bF),max(bF)), (max(bF)-min(bF))/binSize)
    for (i in 2:length(binBorders)) {
      ind<-which(bF>binBorders[i-1] & bF<=binBorders[i])
      keep<-data.frame(c("Chr"=as.character(thisChr[max(ind),"Chr"]),
              thisChr[max(ind),c("Pos","Position")],
              colSums(thisChr[ind,countCols])))
      if (length(newTable)!=0) {
        newTable<-rbind(newTable,keep)
      } else {
        newTable<-keep
      }
    }
  }
  row.names(newTable)<-NULL
  return(newTable)
}



################## functions to evaluate peaks ##############
getChr<-function(cumulativePosition,chrEnds,chrNames) {
  ind<-max(which(cumulativePosition>chrEnds))
  return(chrNames[ind])
}

getChrWidth<-function(cumulativePosition,chrEnds) {
  ind<-max(which(cumulativePosition>chrEnds))
  return(diff(chrEnds)[ind])
}

getValAtLocus<-function(myData, locusPosition, yvals="Gprime") {
  posBefore<-max(which(myData$Position<min(locusPosition)))
  posAfter<-min(which(myData$Position>max(locusPosition)))
  avrSignal<-mean(myData[c(posBefore,posAfter),yvals])
  return(avrSignal)
}

getPeakWidth<-function(myData, myPeak, yvals="Gprime") {
  #calculate peak width at half max
  halfMax<-myPeak$y/2
  ltHalfMax<-which(myData[,yvals]<halfMax)
  peakPos<-myPeak$x
  lowerLim<-max(ltHalfMax[myData$Position[ltHalfMax]<peakPos],1,na.rm=TRUE)
  upperLim<-min(ltHalfMax[myData$Position[ltHalfMax]>peakPos],dim(myData)[1],na.rm=TRUE)
  peakW<-data.frame(pW=upperLim-lowerLim,lhs=lowerLim,rhs=upperLim)
  return(peakW)
}

getPeakQuality<-function(myData,myPeak,locusPosition,chrEnds,yvals="Gprime") {
  peakPos<-myPeak$x
  chrW<-getChrWidth(myPeak$x,chrEnds)
  peakW<-getPeakWidth(myData,myPeak,yvals) #full width at half max
  locH<-getValAtLocus(myData,locusPosition,yvals)
  peakH<-myPeak$y
  # calculate deltaH, deltaPos and deltaW and assym
  deltaH<-abs(peakH-locH)/peakH #difference of yvals at max and at locus range: 0,1
  deltaPos<-min(abs(peakPos-locusPosition))/chrW #distance of locus from max range: -chrW,0,1
  deltaW<-peakW$pW/chrW #peak width as a fraction of chromosome width. range: 0,1
  lw<-(locusPosition-peakW$lhs)
  rw<-(peakW$rhs-locusPosition)
  assym<-ifelse(lw<rw,lw/rw,rw/lw)
  # calculate quality score:
  quality1<-(1-deltaH)*(1-deltaPos)*(1-(myPeak$w/chrW)) #quality score using peaks function width
  quality2<-(1-deltaH)*(1-deltaPos)*(1-deltaW) #quality score using my function for width
  quality3<-(1-deltaH)*(1-deltaPos)*(1-deltaW)*assym #includes penalty for assymetric peaks
  return(c(quality1,quality2,quality3))
}

rr101byChr<-function(myData,recRate,halfWinSize,...) {
  library(DescTools,quietly=TRUE)
  chrs<-levels(myData$Chr)
  #countCols<-grep("counts",names(myData))
  rr101<-c()
  #allGvals<-c()
  for (c in chrs) {
    position<-1:(length(myData$Position[myData$Chr==c])+2*halfWinSize)
    rr<-c(rep(0,halfWinSize),recRate[myData$Chr==c],rep(0,halfWinSize))
    for (i in (halfWinSize+1):(length(position)-halfWinSize)) {
      rr101<-c(rr101,mean(rr[(i-halfWinSize):(i+halfWinSize)]))
    }
  }
  return(rr101)
}
