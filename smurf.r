args <- commandArgs(trailingOnly = TRUE)
resdir <- as.character(args[1])
setwd(resdir)

method<-as.character(args[2])

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

#peaks file containing the coordinates of the regions of interest and their annotation to a gene if promoter (based on gencode) or as DistalRE
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

#calculate global background mutation rate
bmr.global<-sum(mutatedsamples.counts$tot.mut.inRegions)/sum(peaks$lengths) 

#calculate sample-specific background mutation rate
mutatedsamples.counts<-transform(mutatedsamples.counts,bmr.sample=tot.mut.inRegions/sum(peaks$lengths))


#####PERFORM BINOMIAL TEST#####

# x	number of successes, or a vector of length 2 giving the numbers of successes and failures, respectively = FreqUnique
# n	number of trials =peakLEN
# p	hypothesized probability of success =bmr.global or bmr.peak
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
  bmr.peak<-mean(mutatedsamples.counts[which(mutatedsamples.counts[,1] %in% s.in.peak),5]) #average of the mutation rates of the samples with a mutation in that peak
  myd$MutSamples[i]<-as.character(paste(s.in.peak,collapse=",")) 
  myd$MutSamplesunique[i]<-as.character(paste(unique(s.in.peak),collapse=",")) #the list of unique MutSamples samples with muts in this region
  myd$chr[i]<-df[which(df[,17]==myd[i,1])[1],1]
  myd$start[i]<-df[which(df[,17]==myd[i,1])[1],2]
  myd$end[i]<-df[which(df[,17]==myd[i,1])[1],3]
  myd$FreqUnique[i]<-length(unique(s.in.peak)) #the number of samples in that list of unique samples
  myd$peakLEN[i]<-peaks[which(peaks[,7]==myd[i,1]),6]
  myd$peakANNO[i]<-peaks[which(peaks[,7]==myd[i,1]),5]
  myd$peakMUTr[i]=myd[i,2]/myd[i,4]
  if (method=="global") {myd$pbin[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=bmr.global,alternative="greater")$p.value}
  if (method=="peak") {myd$pbin[i]<-binom.test(x=myd$FreqUnique[i],n=myd$peakLEN[i],p=bmr.peak,alternative="greater")$p.value}	  
}



myd<-transform(myd,pbin.adj=p.adjust(myd$pbin,method="BH",n=nrow(peaks))) #use the full number of regions as n for pvalue correction 
myd<-transform(myd,neglog10qval=-log10(myd$pbin.adj))
myd<-myd[order(myd$neglog10qval,decreasing=T),]

myd3<-myd[which(myd$FreqUnique>=3),] #keep only regions with 3 or more mutations from 3 or more samples
myd0.05<-myd3[which(myd3$pbin.adj<=0.05),]


#write to file
write.table(myd,"Mutated_Regions.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd3,"Mutated_Regions_Freq3.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd0.05,"Mutated_Regions_Freq3_qval0.05.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(myd[,6:8],"Mutated_Regions_allmut_DistalREAndProm.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(myd0.05[,6:8],"Mutated_Regions_Freq3_qval0.05_DistalREAndProm.bed",sep="\t",col.names=F,row.names=F,quote=F)

#Make QQplot
pdf("QQplot_Freq3.pdf",height=6,width=9)
qqnorm(-log10(myd3$pbin),main="QQplot",pch=16)
qqline(-log10(myd3$pbin),col="red")
dev.off()

#Write to file some summary stats:
mysum<-paste('myRegions:\t',nrow(peaks),'\n',
             'Regions annotated to Distal REs:\t',nrow(peaks[which(peaks$annotation == "DistalRE"),]),'\n',
             'Regions annotated to Promoters:\t',nrow(peaks[which(peaks$annotation != "DistalRE"),]),'\n',
             'Average number of mutations per sample:\t',mean(mutatedsamples.counts$tot.mut),'\n',
             'Average percentage of mutations in regions:\t',mean(mutatedsamples.counts$perc.mut.inRegions),'\n',
             'Mutated Regions:\t',nrow(myd),'\n',
             'Mutated Regions annotated to Distal REs:\t',nrow(myd[which(myd$peakANNO == "DistalRE"),]),'\n',
             'Mutated Regions annotated to Promoters:\t',nrow(myd[which(myd$peakANNO != "DistalRE"),]),'\n',
             'Mutated Regions (Freq>=3):\t',nrow(myd3),'\n',
             'Mutated Regions (Freq>=3,qval<=0.05):\t',nrow(myd3[which(myd3$pbin.adj<=0.05),]),'\n',
             'Mutated Regions (Freq>=3,qval<=0.05) annotated to Distal REs:\t',nrow(myd3[which(myd3$pbin.adj<=0.05 & myd3$peakANNO == "DistalRE"),]),'\n',
             'Mutated Regions (Freq>=3,qval<=0.05) annotated to Promoters:\t',nrow(myd3[which(myd3$pbin.adj<=0.05 & myd3$peakANNO != "DistalRE"),]),'\n',
             sep="")
write.table(mysum,file="summary.txt",quote=F,row.names=F,col.names=F)

#generate bed file of sig mutated regions annotated to DistalRE for input into c3d
forc3d<-myd0.05[which(myd0.05$peakANNO == "DistalRE"),6:8]
write.table(forc3d,"Mutated_Regions_Freq3_qval0.05_DistalRE.bed",sep="\t",col.names=F,row.names=F,quote=F)


#plot in black and white the peaks pvalue vs mutation rate
if (nrow(myd0.05)>0) {
	pdf("Mutated_Regions_Freq3_qval0.05.pdf",height=8,width=12)
	mymax<-ceiling(max(myd0.05$neglog10qval)/10)*10
	plot(myd0.05$peakMUTr,myd0.05$neglog10qval,ylab="-log10(Qvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax),pch=16,lwd=3)
	dev.off()
}

#Plot
if (nrow(myd0.05)>0) {
	myd0.052<-transform(myd0.05, colour = rep("palevioletred",nrow(myd0.05)),stringsAsFactors=FALSE)
	myd0.052$colour[which(myd0.052$peakANNO == "DistalRE")]<-"lightblue3"

	pdf("Mutated_Regions_Freq3_qval0.05_Promoters_VS_DistalREs.pdf",height=8,width=12)
	numEnh<-nrow(myd0.052[which(myd0.052$peakANNO == "DistalRE"),])
	numProm<-nrow(myd0.052[which(myd0.052$peakANNO != "DistalRE"),])
	mymax<-ceiling(max(myd0.052$neglog10qval)/10)*10
	plot(myd0.052$peakMUTr,myd0.052$neglog10qval,ylab="-log10(Qvalue)",xlab="Region Mutation Rate",ylim=c(0,mymax), col=myd0.052$colour,pch=16,lwd=3,font.lab=2)
	legend(x="topleft", legend = c(paste("Promoter (",numProm,")",sep=""),paste("Distal RE (",numEnh,")",sep="")), col=c("palevioletred","lightblue3"), pch=19, 	cex=1.2, bty="n", pt.cex=1.2,y.intersp=1.4)
	dev.off()
}
