# script to process counts from megamate sequencing. Requires functions in MegaMatePlot.R
# modified 20160129 to make it generalizable for non RNAi phenotypes
# added modifications of different genomes
# modified 20160629 to add latest smoothing functions
#modified 20160630 to use separate G and smG functions

library(signal)
library(RColorBrewer)
library(zoo)
library(dplyr)
options(scipen=10)

source("./MegaMatePlot.R")
genomeVer="WS230"
# fileNames<-list.files(paste0("../finalData/",genomeVer),pattern="vars.*txt")
# samples<-read.table("sampleList.txt",header=TRUE,stringsAsFactors = FALSE)
# codes<-as.vector(t(data.frame(strsplit(fileNames,"_"))[2,]))
# fileNames<-data.frame(fileName=fileNames,code=codes,stringsAsFactors = FALSE)
# fileList<-left_join(fileNames,samples,by="code")
# fileList<-arrange(fileList,phenotype,population,treatment)
# i<-which(names(fileList)!="code")
# fileList<-fileList[,c(names(fileList)[i],"code")]
# write.table(fileList,paste0("../finalData/",genomeVer,"/fileList.txt"),quote=FALSE)

# read in list of experiments and reference genome index
fileList<-read.table(paste0("../finalData/",genomeVer,"/fileList_20160630.txt"),header=TRUE, colClasses="character")
genomeIdx<-read.table(file=paste0("../genomeSeq/",genomeVer,"/c_elegans.",genomeVer,".genomic.fa.fai"))

#load granges object with genomic position of genes in this genome build
load(file=paste0("../../GenomeBuilds/",genomeVer,"/Granges_",genomeVer,".RData"))

# phenotypes<-list(RNAi_mel26=c("zeel-1","ppw-1","mel-26","fem-1"),
#                  RNAi_his55=c("zeel-1","ppw-1","fem-1","his-55"),
#                  Longevity=c("zeel-1","npr-1","glb-5","exp-1","nath-10","rec-8","daf-2","daf-16","age-1","akt-1","daf-18","cki-1","unc-31"),
#                  L1_survival=c("zeel-1","npr-1","glb-5","exp-1","nath-10","lin-4","daf-2","daf-16","age-1","akt-1","daf-18","cki-1","unc-31"))

phenotypes<-list(RNAi_mel26=c("zeel-1","ppw-1","mel-26","fem-1"),
                 RNAi_his55=c("zeel-1","ppw-1","fem-1","his-55"),
                 Longevity=c("zeel-1","npr-1","glb-5","exp-1","nath-10","rec-8"),
                 L1_survival=c("zeel-1","npr-1","glb-5","exp-1","nath-10"))


# bundle repatively used data into list
chrs=c("CHROMOSOME_I","CHROMOSOME_II","CHROMOSOME_III","CHROMOSOME_IV","CHROMOSOME_V","CHROMOSOME_X")
chrNum<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX")
#locNames<-c()#locNames<-as.character(gff_a[gff_a$public %in% phenotypes$Longevity,]$public)
#locs=c()#locs<-c(2342232,4185076,9132545,50551562,56356262)
labelAdj<-c(0)#rep(0,length(locs))#c(0.9,0.2,0.1,NA,NA)
chrEnds<-c(0,cumsum(genomeIdx$V2))[1:(length(genomeIdx$V2))]
midpoints<-chrEnds[1:(length(chrEnds)-1)]+diff(chrEnds)/2
bundle<-list(chrs=chrs, chrNum=chrNum, locNames=c(), locs=c(), labelAdj=labelAdj, chrEnds=chrEnds,
             midpoints=midpoints, mainTitle="",legend="")

### to find distance between two gr objects (my peaks and the possible loci)
# locOfInterest<-gff_a[gff_a$public %in% phenotypes$RNAi_mel26,]
# peak=GRanges(seqnames=c("CHROMOSOME_I"),
#            ranges=IRanges(start=c(5000000),end=c(5001000)),
#            strand=c("+"))
# locdistance<-distanceToNearest(peak,locOfInterest)
# mcols(locdistance)$distance
# locOfInterest[subjectHits(locdistance),]
#########

d<-format(Sys.time(), "%Y%m%d")

# #create a column in which to write name of file with processed data
# fileList["processedData"]<-NA

for (ff in seq(1,dim(fileList)[1],by=2)) { 
  # # read in control and selected data for single experiment
  # contFile=fileList[fileList$treatment=="control",1][(ff+1)/2]
  # selectFile=fileList[fileList$treatment=="selected",1][(ff+1)/2]
  # cont<-read.table(paste0("../finalData/",genomeVer,"/",contFile),header=TRUE)
  # select<-read.table(paste0("../finalData/",genomeVer,"/",selectFile),header=TRUE)

  # read in processed data
  myExpFile=fileList[fileList$treatment=="control","processedData"][(ff+1)/2]
  myExp<-read.table(myExpFile,header=TRUE)
  
  # # add cumulative genome coordinates to data table
  # cont<-cumPosition(cont,genomeIdx)
  # select<-cumPosition(select,genomeIdx)

  # # clean data by removing data with less than minRead, and whose read numbers are extreme outliers
  # cont<-cleanData(cont, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]), genomeVer=genomeVer,
  #                 MADs=2, minReads=5, DelNonPolymorphic=FALSE)
  # select<-cleanData(select, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]), genomeVer=genomeVer,
  #                   MADs=2, minReads=5, DelNonPolymorphic=FALSE)
  # # merge the two tables on common columns (chr, position, variant)
  # myExp<-merge(cont,select)
  # # remove loci that are not polymorphic in both control and selected
  # myExp<-removeNonPolymorphic(myExp, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]),genomeVer=genomeVer)
  # # remove superfluous columns (separate counts of forward and reverse reads)
  # excl<-grep("fwd|rev",names(myExp))
  # myExp<-myExp[,-excl]
  # # order table according to genome coordinates
  # myExp<-myExp[order(myExp$Position),]
  
  # # do fisher.test for each locus
  # countCols<-grep("counts",names(myExp))
  # myExp["fisherExact.pvals"]<-sapply(1:dim(myExp)[1],function(i){
  #   fisher.test(matrix(as.matrix(myExp[i,countCols]),ncol=2))$p.value})
  # 
  # # adjust p values for multiple testing
  # myExp["padj.fdr"]<-p.adjust(myExp[,"fisherExact.pvals"],method="fdr")
  # 
  # # do GTest
  # G<-getGbyChr(myExp,myExp$Position)
  # myExp<-cbind(myExp,Gval.GTest=G)
  #smG<-smootheGbyChr(myExp,rep(1000000,dim(myExp)[1]),myExp$Position,"Gval.GTest")
  #myExp<-cbind(myExp,Gval.GTest=G)
  
  #do GTest with python function
  #First write a file with counts
  #write.table
  
  ########
  # open pdf for plots
  pdf(file=paste0("../finalData/",genomeVer,"/",d,"_",fileList[ff,"population"],".pdf"), paper="a4", height=11, width=8)
  par(mfrow=c(4,1))
  bundle$mainTitle<-substitute(paste(x," ",y), list(x=fileList[ff,"population"],y=fileList[ff,"phenotype"]))
  bundle$locNames<-getLocNames(gff_a,phenotypes,fileList[ff,"phenotype"])
  bundle$locs<-getLocPositions(gff_a,phenotypes,fileList[ff,"phenotype"],genomeIdx)
  useLoci<-1:length(bundle$locs)

  # plot raw Hawaii allele frequencies of control and selected along genome
  plotFreq(myExp,extraData=bundle,useLoci=useLoci)
  
  # plot the difference between Hawaii allele frequencies in control and selected populations
  # smoothe with sgolay polynomial smoothing (window size=sm) before plotting
  plotSmoothedDiff(myExp,extraData=bundle,sm=1001,useLoci=useLoci)
    
  # plot -log10 of p values after polynomial smoothing (sgolay filter)
  plotSmPvals(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
  
  #plot sgolay smoothed G
  plotSmGvals_sGolay(myExp, myExp$Gval.GTest, extraData=bundle, sm=1001, useLoci=useLoci)
  
  # #plot smG
  # chrList<-levels(myExp$Chr)
  # plot(myExp$Position,smG,type='n',lwd=2,ylab="Smoothed G value",
  #      xlab="Position (bp)",main=bundle$mainTitle,ylim=c(0,max(smG)*1.1))
  # plotByChr(myExp, smG, chrList, chrColumn="Chr",lwd=2)
  # 
  # text(bundle$midpoint,max(smG)*1.1,bundle$chrNum,cex=1, col="black")
  # abline(v=bundle$chrEnds,col="dark gray",lty=5)
  # abline(v=bundle$locs[useLoci], col="red",lwd=0.8)
  # mtext(bundle$locNames[useLoci], at=bundle$locs[useLoci], cex=0.7, col="red",
  #       adj=bundle$labelAdj[useLoci])
  # title(sub="2Mb window smoothed GTest statistic (1000 permutations)")
  dev.off()
  ##########
    
  # #output processed data as a table
  # processedFileName<-paste0("../finalData/",genomeVer,"/processedData_",d,"_",fileList[ff,"phenotype"],
  #                          "_",fileList[ff,"population"],".txt")
  # write.table(myExp,file=processedFileName)
  
  # #add filename of processedData to fileList table
  # fileList[ff,"processedData"]<-processedFileName
  # fileList[ff+1,"processedData"]<-processedFileName
}

# #write new fileList_date.txt file to store processedFileName 
# write.table(fileList,file=paste0("../finalData/",genomeVer,"/fileList_",d,".txt"))

plot=FALSE
if (plot==TRUE) {
#plot histograms of readDepths
pdf(file=paste0("../finalData/",genomeVer,"/",d,"_readDepths.pdf"), paper="a4", height=11, width=8)
par(mfrow=c(4,2))
for (ff in seq(1,dim(fileList)[1],by=2)) { 
  # read in control and selected data for single experiment
  myExp<-read.table(as.character(fileList[ff,"processedData"]),header=TRUE)
  
  # plot histograms of counts
  rD<-grep("readDepth",names(myExp))
  hist(myExp[,rD[1]],breaks=100, main=paste("Control",fileList[ff,2],fileList[ff,3]),xlab="Read depth")
  abline(v=mean(myExp[,rD[1]]),col="red")
  hist(myExp[,rD[2]],breaks=100, main=paste("Selected",fileList[ff,2],fileList[ff,3]),xlab="Read depth")
  abline(v=mean(myExp[,rD[2]]),col="red")
  
  #remove nonPolymorphic loci
  myExp<-removeNonPolymorphic(myExp, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]),genomeVer=genomeVer)
  
  # plot histograms of readDepths
  rD<-grep("readDepth",names(myExp))
  hist(myExp[,rD[1]],breaks=100, main=paste("Control",fileList[ff,2],fileList[ff,3]),xlab="Read depth")
  abline(v=mean(myExp[,rD[1]]),col="red")
  hist(myExp[,rD[2]],breaks=100, main=paste("Selected",fileList[ff,2],fileList[ff,3]),xlab="Read depth")
  abline(v=mean(myExp[,rD[2]]),col="red")
}
dev.off()
}

