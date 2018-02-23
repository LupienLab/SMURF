
# This file is part of SMuRF
# Copyright 2017 Paul Guilhamon

# SMuRF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# SMuRF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details <http://www.gnu.org/licenses/>.



#######################################################
# Written by: Paul Guilhamon
# Princess Margaret Cancer Centre - University Health Network, November 2017
#######################################################







args <- commandArgs(trailingOnly = TRUE)
resdir <- as.character(args[1])
setwd(resdir)

method<-as.character(args[2])


####LOAD IN FILES & PACKAGES####
suppressMessages(library(gtools))
suppressMessages(library(gplots))
suppressMessages(library(data.table))
suppressMessages(library(psych))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))


`%ni%` <- Negate(`%in%`) 

#full VCF file containing all somatic variants
df<-read.delim("FiltVarsInPeaks.txt",sep="\t",header=F,stringsAsFactors = FALSE)
names(df)<-c("chr","start","end","samples","MutSamples","chr_mut","pos_mut")

#peaks file containing the coordinates of the regions of interest and their annotation to a gene if promoter (based on gencode) or as DistalRE
peaks<-read.delim("annotated_regions.txt",sep="\t",header=F,stringsAsFactors=FALSE)
names(peaks)<-c("chr","start","end","samples","annotation")

#number of SNV in each sample
mutatedsamples.counts<-read.delim("named_counts.txt",sep="\t",header=F,stringsAsFactors=FALSE)
names(mutatedsamples.counts)<-c("MutSamples","tot.mut")


#from Rose code (https://bitbucket.org/young_computation/rose/src/1a9bb86b546476d0e227640aa2995332a3b650c3/ROSE_callSuper.R?at=master&fileviewer=file-view-default):
#X11 License
#Copyright (C) 2013 by Whitehead Institute for Biomedical Research

#This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
  inputVector <- sort(inputVector)
  #inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  if(drawPlot){  #if TRUE, draw the plot
    plot(1:length(inputVector), inputVector,type="l",...)
    b <- y_cutoff-(slope* xPt)
    abline(v= xPt,h= y_cutoff,lty=2,col=8)
    points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
    abline(coef=c(b,slope),col=2)
    title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
  }
  return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}
#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}


####CALCULATE REQUIRED INFO FOR BINOMIAL TEST####

#add in number of muts per sample over the genomic regions
mutatedsamples.counts<-transform(mutatedsamples.counts,tot.mut.inRegions=rep(0,nrow(mutatedsamples.counts)))
for (i in 1:nrow(mutatedsamples.counts)){
  mutatedsamples.counts[i,3]<-nrow(df[which(df$MutSamples==mutatedsamples.counts[i,1]),])
}

mutatedsamples.counts<-transform(mutatedsamples.counts,perc.mut.inRegions=100*(tot.mut.inRegions/tot.mut))

pdf("barplot_mut_inRegions_Counts_and_Perc.pdf",height=6,width=9)
par(mar=c(8, 4, 4, 2))
mutatedsamples.counts<-mutatedsamples.counts[order(mutatedsamples.counts$tot.mut,decreasing=T),]
mymax1<-ceiling(max(mutatedsamples.counts$tot.mut)/5000)*5000
mymax2<-ceiling(max(mutatedsamples.counts$perc.mut.inRegions)/5)*5
barplot(mutatedsamples.counts$tot.mut,names.arg=mutatedsamples.counts$MutSamples,las=2,ylim=c(0,mymax1),cex.names=0.8,ylab="Total number SNVs",col="#fdcdac")
abline(h=mean(mutatedsamples.counts$tot.mut),col="red")
barplot(mutatedsamples.counts$perc.mut.inRegions,names.arg=mutatedsamples.counts$MutSamples,las=2,ylim=c(0,mymax2),cex.names=0.8,ylab="%SNVs in myRegions",col="#b3e2cd")
abline(h=mean(mutatedsamples.counts$perc.mut.inRegions),col="red")
dev.off()


#calculate length of each peak
peaks<-transform(peaks, lengths=end-start)

#calculate number of mutations in each region
peaks<-transform(peaks,peakID=paste(paste(as.character(chr),as.character(start),sep=":"),as.character(end),sep="-"))
df<-transform(df,peakID=paste(paste(as.character(chr),as.character(start),sep=":"),as.character(end),sep="-"))

peaks.counts<-as.data.frame(table(df$peakID),stringsAsFactors=F) 
names(peaks.counts)<-c("peakID","Freq")
peaks.counts<-peaks.counts[order(peaks.counts$Freq,decreasing=T),]


#calculate sample-specific background mutation rate
PeaksWithMuts<-peaks[which(peaks$peakID %in% unique(df$peakID)),]
mutatedsamples.counts<-transform(mutatedsamples.counts,bmr.sample=tot.mut.inRegions/sum(PeaksWithMuts$lengths))

#calculate allsamples background mutation rate
bmr.allsamples<-mean(mutatedsamples.counts[,"bmr.sample"])


#####PERFORM BINOMIAL TEST#####

# x	number of successes, or a vector of length 2 giving the numbers of successes and failures, respectively = FreqUnique
# n	number of trials =peakLEN
# p	hypothesized probability of success =bmr.allsamples or bmr.regionsamples
# alternative	indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter.: "greater"
# conf.level confidence level for the returned confidence interval:0.95

myd<-transform(peaks.counts,
               FreqUnique=rep(0,nrow(peaks.counts)),
               peakLEN=rep(0,nrow(peaks.counts)),
               peakANNO=rep("bla",nrow(peaks.counts)),
               chr=rep("chr",nrow(peaks.counts)),start=rep(0,nrow(peaks.counts)),end=rep(0,nrow(peaks.counts)),
               MutSamples=rep("bla",nrow(peaks.counts)),
               MutSamplesunique=rep("bla",nrow(peaks.counts)),
               peakMUTr=rep(0,nrow(peaks.counts)),
               pbin_as=rep(0,nrow(peaks.counts)),
               pbin_rs=rep(0,nrow(peaks.counts)),stringsAsFactors = FALSE)


for (i in 1:nrow(myd)){
  s.in.peak<-df[which(df[,"peakID"]==myd[i,"peakID"]),"MutSamples"] #which MutSamples samples have mutations in that region
  bmr.regionsamples<-mean(mutatedsamples.counts[which(mutatedsamples.counts[,"MutSamples"] %in% s.in.peak),"bmr.sample"]) #average of the mutation rates of the samples with a mutation in that peak
  myd$FreqUnique[i]<-length(unique(s.in.peak)) #the number of unique mutated samples in the peak
  myd$peakLEN[i]<-peaks[which(peaks[,"peakID"]==myd[i,"peakID"]),"lengths"]
  myd$peakANNO[i]<-peaks[which(peaks[,"peakID"]==myd[i,"peakID"]),"annotation"]
  myd$chr[i]<-df[which(df[,"peakID"]==myd[i,"peakID"])[1],"chr"]
  myd$start[i]<-df[which(df[,"peakID"]==myd[i,"peakID"])[1],"start"]
  myd$end[i]<-df[which(df[,"peakID"]==myd[i,"peakID"])[1],"end"]
  myd$MutSamples[i]<-as.character(paste(s.in.peak,collapse=",")) 
  myd$MutSamplesunique[i]<-as.character(paste(unique(s.in.peak),collapse=",")) #the list of unique MutSamples samples with muts in this region
  myd$peakMUTr[i]<-myd[i,"FreqUnique"]/myd[i,"peakLEN"]
  #if (method=="allsamples") {myd$pbin[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=bmr.allsamples,alternative="greater")$p.value}
  #if (method=="regionsamples") {myd$pbin[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=bmr.regionsamples,alternative="greater")$p.value}
  myd$pbin_as[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=bmr.allsamples,alternative="greater")$p.value
  myd$pbin_rs[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=bmr.regionsamples,alternative="greater")$p.value
}




myd<-transform(myd,pbin_as.adj=p.adjust(myd$pbin_as,method="BH",n=nrow(peaks)),stringsAsFactors = FALSE) #use the full number of regions as n for pvalue correction
myd<-transform(myd,pbin_rs.adj=p.adjust(myd$pbin_rs,method="BH",n=nrow(peaks)),stringsAsFactors = FALSE) #use the full number of regions as n for pvalue correction 
myd<-transform(myd,neglog10qval_as=-log10(myd$pbin_as.adj),stringsAsFactors = FALSE)
myd<-transform(myd,neglog10qval_rs=-log10(myd$pbin_rs.adj),stringsAsFactors = FALSE)

if (method=="allsamples") {
  myd<-myd[order(myd$neglog10qval_as,decreasing=T),]
  mydsig<-myd[which(myd$pbin_as.adj<=0.05),]
  }

if (method=="regionsamples") {
  myd<-myd[order(myd$neglog10qval_rs,decreasing=T),]
  mydsig<-myd[which(myd$pbin_rs.adj<=0.05),]
  }

mythreshold<-calculate_cutoff(mydsig$peakMUTr,drawPlot=F)$absolute
mydFINAL<-mydsig[which(mydsig$peakMUTr>=mythreshold),]




#write to file
write.table(myd,"Mutated_Regions.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(mydsig,"Mutated_Regions_qval0.05.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(mydFINAL,file=paste0("Mutated_Regions_qval0.05_peakMUTr",mythreshold,".txt"),sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd[,6:8],"Mutated_Regions_allmut_DistalREAndProm.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(mydFINAL[,6:8],file=paste0("Mutated_Regions_qval0.05_peakMUTr",mythreshold,"_DistalREAndProm.bed"),sep="\t",col.names=F,row.names=F,quote=F)



#Write to file some summary stats:
mysum<-paste('myRegions:\t',nrow(peaks),'\n',
             'Regions annotated to Distal REs:\t',nrow(peaks[which(peaks$annotation == "DistalRE"),]),'\n',
             'Regions annotated to Promoters:\t',nrow(peaks[which(peaks$annotation != "DistalRE"),]),'\n',
             'Average number of mutations per sample:\t',mean(mutatedsamples.counts$tot.mut),'\n',
             'Average percentage of mutations in regions:\t',mean(mutatedsamples.counts$perc.mut.inRegions),'\n',
             'Mutated Regions:\t',nrow(myd),'\n',
             'Mutated Regions annotated to Distal REs:\t',nrow(myd[which(myd$peakANNO == "DistalRE"),]),'\n',
             'Mutated Regions annotated to Promoters:\t',nrow(myd[which(myd$peakANNO != "DistalRE"),]),'\n',
             'Mutated Regions (qval<=0.05):\t',nrow(mydsig),'\n',
             'Mutated Regions (qval<=0.05, peakMUTr>=threshold):\t',nrow(mydFINAL),'\n',
             'Mutated Regions (qval<=0.05, peakMUTr>=threshold) annotated to Distal REs:\t',nrow(mydFINAL[which(mydFINAL$peakANNO == "DistalRE"),]),'\n',
             'Mutated Regions (qval<=0.05, peakMUTr>=threshold) annotated to Promoters:\t',nrow(mydFINAL[which(mydFINAL$peakANNO != "DistalRE"),]),'\n',
             sep="")
write.table(mysum,file="summary.txt",quote=F,row.names=F,col.names=F)

#generate bed file of sig mutated regions annotated to DistalRE for input into c3d
forc3d<-mydFINAL[which(mydFINAL$peakANNO == "DistalRE"),6:8]
write.table(forc3d,file=paste0("Mutated_Regions_qval0.05_peakMUTr",mythreshold,"_DistalRE.bed"),sep="\t",col.names=F,row.names=F,quote=F)



#####PLOT ADJUSTED P-VALUE vs REGION MUTATION RATE#####
#plot in black and white the peaks pvalue vs mutation rate
if (nrow(mydFINAL)>0) {
  pdf(file=paste0("Mutation_Rate_Plot_qval0.05_peakMUTr",mythreshold,".pdf"),height=8,width=12)
  if (method=="allsamples") {
    mymax<-ceiling(max(mydFINAL$neglog10qval_as)/10)*10
    plot(mydFINAL$peakMUTr,mydFINAL$neglog10qval_as,ylab="-log10(Qvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax),pch=16,lwd=3)
  }
  else if (method=="regionsamples") {
    mymax<-ceiling(max(mydFINAL$neglog10qval_rs)/10)*10
    plot(mydFINAL$peakMUTr,mydFINAL$neglog10qval_rs,ylab="-log10(Qvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax),pch=16,lwd=3)
  }
  dev.off()
}

#Plot in colour, separating promoters and Distal REs
if (nrow(mydFINAL)>0) {
  mydFINALcol<-transform(mydFINAL, colour = rep("palevioletred",nrow(mydFINAL)),stringsAsFactors=FALSE)
  mydFINALcol$colour[which(mydFINALcol$peakANNO == "DistalRE")]<-"lightblue3"
  
  pdf(file=paste0("Mutation_Rate_Plot_qval0.05_peakMUTr",mythreshold,"_Promoters_VS_DistalREs.pdf"),height=8,width=12)
  numEnh<-nrow(mydFINALcol[which(mydFINALcol$peakANNO == "DistalRE"),])
  numProm<-nrow(mydFINALcol[which(mydFINALcol$peakANNO != "DistalRE"),])
  if (method=="allsamples") {
    mymax<-ceiling(max(mydFINALcol$neglog10qval_as)/10)*10
    plot(mydFINALcol$peakMUTr,mydFINALcol$neglog10qval_as,ylab="-log10(Qvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax), col=mydFINALcol$colour,pch=16,lwd=3,font.lab=2)
  }
  else if (method=="regionsamples") {
    mymax<-ceiling(max(mydFINALcol$neglog10qval_rs)/10)*10
    plot(mydFINALcol$peakMUTr,mydFINALcol$neglog10qval_rs,ylab="-log10(Qvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax), col=mydFINALcol$colour,pch=16,lwd=3,font.lab=2)
  }
  legend(x="topleft", legend = c(paste("Promoter (",numProm,")",sep=""),paste("Distal RE (",numEnh,")",sep="")), col=c("palevioletred","lightblue3"), pch=19, 	cex=1.2, bty="n", pt.cex=1.2,y.intersp=1.4)
  dev.off()
}








#####PLOT ADJUSTED P-VALUE vs FREQUENCY OF UNIQUE SAMPLES MUTATED IN REGION#####

#in black and white qvalue vs Mutation frequency
if (method=="allsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_as)) + 
    geom_point(aes(color = "black"), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="none") + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = "#000000") # change colors
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,".pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_as)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  bwplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(bwplot)
  dev.off()
}

if (method=="regionsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_rs)) + 
    geom_point(aes(color = "black"), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="none") + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = "#000000") # change colors
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,".pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_rs)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  bwplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(bwplot)
  dev.off()
}





#in colours qvalue vs Mutation frequency, with distalREs and Promoters
numEnh<-nrow(mydFINAL[which(mydFINAL$peakANNO=="DistalRE"),])
numProm<-nrow(mydFINAL[which(mydFINAL$peakANNO!="DistalRE"),])


if (method=="allsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL %>%
    mutate(color = ifelse(mydFINAL$peakANNO!="DistalRE",
                          yes = "Promoter",
                          no = "DistalRE"))
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_as)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="bottom",legend.title = element_blank(),legend.key = element_rect(colour = NA)) + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = c("DistalRE"="lightblue3","Promoter"="palevioletred"),
                       labels=c("DistalRE"=paste0("Distal RE (",numEnh,")"),"Promoter"=paste0("Promoter (",numProm,")")))
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,"_Promoters_VS_DistalREs.pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_as)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  colplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(colplot)
  dev.off()
}

if (method=="regionsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL %>%
    mutate(color = ifelse(mydFINAL$peakANNO!="DistalRE",
                          yes = "Promoter",
                          no = "DistalRE"))
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_rs)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="bottom",legend.title = element_blank(),legend.key = element_rect(colour = NA)) + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = c("DistalRE"="lightblue3","Promoter"="palevioletred"),
                       labels=c("DistalRE"=paste0("Distal RE (",numEnh,")"),"Promoter"=paste0("Promoter (",numProm,")")))
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,"_Promoters_VS_DistalREs.pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_rs)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  colplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(colplot)
  dev.off()
}





#Plot of promoters and DistalREs separately
##DistalREs
numEnh<-nrow(mydFINAL[which(mydFINAL$peakANNO=="DistalRE"),])
numProm<-nrow(mydFINAL[which(mydFINAL$peakANNO!="DistalRE"),])



if (method=="allsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL %>%
    mutate(color = ifelse(mydFINAL$peakANNO!="DistalRE",
                          yes = "Promoter",
                          no = "DistalRE"))
  data<-data[which(data$peakANNO=="DistalRE"),]
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_as)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="bottom",legend.title = element_blank(),legend.key = element_rect(colour = NA)) + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = c("DistalRE"="lightblue3","Promoter"="palevioletred"),
                       labels=c("DistalRE"=paste0("Distal RE (",numEnh,")"),"Promoter"=paste0("Promoter (",numProm,")")))
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,"_DistalREs.pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_as)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  colplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(colplot)
  dev.off()
}

if (method=="regionsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL %>%
    mutate(color = ifelse(mydFINAL$peakANNO!="DistalRE",
                          yes = "Promoter",
                          no = "DistalRE"))
  data<-data[which(data$peakANNO=="DistalRE"),]
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_rs)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="bottom",legend.title = element_blank(),legend.key = element_rect(colour = NA)) + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = c("DistalRE"="lightblue3","Promoter"="palevioletred"),
                       labels=c("DistalRE"=paste0("Distal RE (",numEnh,")"),"Promoter"=paste0("Promoter (",numProm,")")))
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,"_DistalREs.pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_rs)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  colplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(colplot)
  dev.off()
}







##Promoters
if (method=="allsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL %>%
    mutate(color = ifelse(mydFINAL$peakANNO!="DistalRE",
                          yes = "Promoter",
                          no = "DistalRE"))
  data<-data[which(data$peakANNO!="DistalRE"),]
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_as)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="bottom",legend.title = element_blank(),legend.key = element_rect(colour = NA)) + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = c("DistalRE"="lightblue3","Promoter"="palevioletred"),
                       labels=c("DistalRE"=paste0("Distal RE (",numEnh,")"),"Promoter"=paste0("Promoter (",numProm,")")))
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,"_Promoters.pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_as)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  colplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(colplot)
  dev.off()
}

if (method=="regionsamples" & nrow(mydFINAL)>0) {
  data <- mydFINAL %>%
    mutate(color = ifelse(mydFINAL$peakANNO!="DistalRE",
                          yes = "Promoter",
                          no = "DistalRE"))
  data<-data[which(data$peakANNO!="DistalRE"),]
  colored <- ggplot(data, aes(x = FreqUnique, y = neglog10qval_rs)) + 
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position="bottom",legend.title = element_blank(),legend.key = element_rect(colour = NA)) + # remove legend 
    ggtitle(label = "") +  # add title
    xlab(expression("Samples Mutated in Region")) + # x-axis label
    ylab(expression(-log[10]("q-value"))) + # y-axis label
    scale_color_manual(values = c("DistalRE"="lightblue3","Promoter"="palevioletred"),
                       labels=c("DistalRE"=paste0("Distal RE (",numEnh,")"),"Promoter"=paste0("Promoter (",numProm,")")))
  
  pdf(file=paste0("Sample_Frequency_Plot_qval0.05_peakMUTr",mythreshold,"_Promoters.pdf"),height=8,width=6)
  mymax<-ceiling(max(data$neglog10qval_rs)/10)*10
  mymaxx<-max(data$FreqUnique)
  if (mymaxx<=10) {myinterval<-1}
  else if (mymaxx>10 & mymaxx<=20) {myinterval<-2}
  else if (mymaxx>20 & mymaxx<=40) {myinterval<-5}
  else if (mymaxx>40) {myinterval<-10}
  colplot<-colored + scale_y_continuous(trans = "log1p",breaks=c(seq(0,10,2),seq(0,mymax,10))) + scale_x_continuous(breaks=seq(0,mymaxx,myinterval))
  print(colplot)
  dev.off()
}





