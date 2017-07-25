args <- commandArgs(trailingOnly = TRUE)
resdir <- as.character(args[1])
setwd(resdir)

####LOAD IN FILES & PACKAGES####
library(gtools)
library(gplots)
library(data.table)
library(psych)
library(GenomicRanges)

`%ni%` <- Negate(`%in%`) 

#full VCF file containing all somatic variants
df<-read.delim("fullVCF.txt",sep="\t",header=F,stringsAsFactors = FALSE)
names(df)<-c("chr","start","end","samples","MutSamples","CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

#peaks file containing the coordinates of the regions of interest and their annotation to a gene if promoter (based on gencode) or as enhancer
peaks<-read.delim("annotated_regions.txt",sep="\t",header=F,stringsAsFactors=FALSE)
names(peaks)<-c("chr","start","end","samples","annotation")

#number of SNV in each sample
mutatedsamples.counts<-read.delim("named_counts.txt",sep="\t",header=F,stringsAsFactors=FALSE)
names(mutatedsamples.counts)<-c("MutSamples","tot.mut")

####CALCULATE REQUIRED INFO FOR BINOMIAL TEST####

#add in number of muts per sample over the genomic regions
mutatedsamples.counts<-transform(mutatedsamples.counts,tot.mut.inRegions=rep(0,nrow(mutatedsamples.counts)))
for (i in 1:nrow(mutatedsamples.counts)){
  mutatedsamples.counts[i,3]<-nrow(df[which(df$MutSamples==mutatedsamples.counts[i,1]),])
}

mutatedsamples.counts<-transform(mutatedsamples.counts,perc.mut.inRegions=100*(tot.mut.inRegions/tot.mut))


pdf("barplot_mut_inRegions_Counts_and_Perc.pdf",height=6,width=9)
par(mar=c(8, 4, 4, 2))
mymax1<-ceiling(max(mutatedsamples.counts$tot.mut)/5000)*5000
mymax2<-ceiling(max(mutatedsamples.counts$perc.mut.inRegions)/5)*5
barplot(mutatedsamples.counts$tot.mut,names.arg=mutatedsamples.counts$MutSamples,las=2,ylim=c(0,mymax1),cex.names=0.8,ylab="SNVs",col="blue")
barplot(mutatedsamples.counts$perc.mut.inRegions,names.arg=mutatedsamples.counts$MutSamples,las=2,ylim=c(0,mymax2),cex.names=0.8,ylab="%SNVs in myRegions",col="beige")
boxplot(mutatedsamples.counts$perc.mut.inRegions,las=2,ylim=c(0,mymax2),cex.names=0.8,ylab="%SNVs in myRegions")
dev.off()

#calculate length of each peak
peaks<-transform(peaks, lengths=end-start)

#calculate number of mutations in each region
peaks<-transform(peaks,peakID=paste(paste(as.character(chr),as.character(start),sep=":"),as.character(end),sep="-"))
df<-transform(df,peakID=paste(paste(as.character(chr),as.character(start),sep=":"),as.character(end),sep="-"))

peaks.counts<-as.data.frame(table(df$peakID),stringsAsFactors=F) 
names(peaks.counts)<-c("peakID","Freq")
peaks.counts<-peaks.counts[order(peaks.counts$Freq,decreasing=T),]

#calculate background mutation rate
fullrate<-sum(mutatedsamples.counts$tot.mut.inRegions)/sum(peaks$lengths) 

#####PERFORM BINOMIAL TEST#####

# x	number of successes, or a vector of length 2 giving the numbers of successes and failures, respectively = FreqUnique
# n	number of trials =peakLEN
# p	hypothesized probability of success =fullrate
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
               pbin=rep(0,nrow(peaks.counts)),stringsAsFactors = FALSE)

for (i in 1:nrow(myd)){
  s.in.peak<-df[which(df[,17]==myd[i,1]),5] #which MutSamples samples have mutations in that region
  myd$MutSamples[i]<-as.character(paste(s.in.peak,collapse=",")) 
  myd$MutSamplesunique[i]<-as.character(paste(unique(s.in.peak),collapse=",")) #the list of unique MutSamples samples with muts in this region
  myd$chr[i]<-df[which(df[,17]==myd[i,1])[1],1]
  myd$start[i]<-df[which(df[,17]==myd[i,1])[1],2]
  myd$end[i]<-df[which(df[,17]==myd[i,1])[1],3]
  myd$FreqUnique[i]<-length(unique(s.in.peak)) #the number of samples in that list of unique samples
  myd$peakLEN[i]<-peaks[which(peaks[,7]==myd[i,1]),6]
  myd$peakANNO[i]<-peaks[which(peaks[,7]==myd[i,1]),5]
  myd$peakMUTr[i]=myd[i,2]/myd[i,4]
  myd$pbin[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=fullrate,alternative="greater")$p.value
}

myd<-transform(myd,pbin.adj=p.adjust(myd$pbin,method="BH",n=nrow(peaks))) #use the full number of regions as n for pvalue correction 
myd<-transform(myd,neglog10qval=-log10(myd$pbin.adj))
myd<-myd[order(myd$neglog10qval,decreasing=T),]

myd3<-myd[which(myd$Freq>=3),] #keep only regions with more than 2 mutations
myd0.05<-myd3[which(myd3$pbin.adj<=0.05),]


#write to file
write.table(myd,"Mutated_Regions.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd3,"Mutated_Regions_Freq3.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd0.05,"Mutated_Regions_Freq3_qval0.05.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd[,6:8],"Mutated_Regions_allmut_EnhAndProm.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(myd0.05[,6:8],"Mutated_Regions_sigmut0.05_EnhAndProm.bed",sep="\t",col.names=F,row.names=F,quote=F)

#Make QQplot
pdf("QQplot_Freq3.pdf",height=6,width=9)
qqnorm(-log10(myd3$pbin),main="QQplot",pch=16)
qqline(-log10(myd3$pbin),col="red")
dev.off()


#Write to file some summary stats:
mysum<-paste('myRegions:\t',nrow(peaks),'\n',
             'Regions annotated to Enhancers:\t',nrow(peaks[which(peaks$annotation == "Enhancer"),]),'\n',
             'Regions annotated to Promoters:\t',nrow(peaks[which(peaks$annotation != "Enhancer"),]),'\n',
             'Mutated Regions:\t',nrow(myd),'\n',
             'Mutated Regions annotated to Enhancers:\t',nrow(myd[which(myd$peakANNO == "Enhancer"),]),'\n',
             'Mutated Regions annotated to Promoters:\t',nrow(myd[which(myd$peakANNO != "Enhancer"),]),'\n',
             'Mutated Regions (Freq>=3):\t',nrow(myd3),'\n',
             'Mutated Regions (Freq>=3,qval<=0.05):\t',nrow(myd3[which(myd3$pbin.adj<=0.05),]),'\n',
             'Mutated Regions (Freq>=3,qval<=0.05) annotated to Enhancers:\t',nrow(myd3[which(myd3$pbin.adj<=0.05 & myd3$peakANNO == "Enhancer"),]),'\n',
             'Mutated Regions (Freq>=3,qval<=0.05) annotated to Promoters:\t',nrow(myd3[which(myd3$pbin.adj<=0.05 & myd3$peakANNO != "Enhancer"),]),'\n',
             sep="")
write.table(mysum,file="summary.txt",quote=F,row.names=F,col.names=F)

#generate bed file of sig mutated regions annotated to enhancer for input into c3d
forc3d<-myd[which(myd$pbin.adj<=0.05 & myd$peakANNO == "Enhancer"),6:8]
write.table(forc3d,"Mutated_Regions_sigmut0.05_Enh_forc3d.bed",sep="\t",col.names=F,row.names=F,quote=F)


#plot in black and white the peaks pvalue vs mutation rate
pdf("hotspotsig0.05.pdf",height=8,width=12)
mymax<-ceiling(max(myd0.05$neglog10qval)/10)*10
plot(myd0.05$peakMUTr,myd0.05$neglog10qval,ylab="-log10(Pvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax),pch=16,lwd=3)
dev.off()


#Plot

myd0.052<-transform(myd0.05, colour = rep("palevioletred",nrow(myd0.05)),stringsAsFactors=FALSE)
myd0.052$colour[which(myd0.052$peakANNO == "Enhancer")]<-"lightblue3"

pdf("hotspotsig0.05_PromVEnh.pdf",height=8,width=12)
numEnh<-nrow(myd0.052[which(myd0.052$peakANNO == "Enhancer"),])
numProm<-nrow(myd0.052[which(myd0.052$peakANNO != "Enhancer"),])
mymax<-ceiling(max(myd0.052$neglog10qval)/10)*10
plot(myd0.052$peakMUTr,myd0.052$neglog10qval,ylab="-log10(Pvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax), col=myd0.052$colour,pch=16,lwd=3,font.lab=2)
legend(x="topleft", legend = c(paste("Promoter (",numProm,")",sep=""),paste("Enhancer (",numEnh,")",sep="")), col=c("palevioletred","lightblue3"), pch=19, cex=1.2, bty="n", pt.cex=1.2,y.intersp=1.4)
dev.off()

