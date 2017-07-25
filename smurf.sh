#!/bin/bash
#$ -cwd
#$ -N smurf
#$ -j y
#$ -q light.q
#$ -S /bin/bash

#SMURF v4.0
#Developed by Paul Guilhamon

#Brief description of arguments; please refer to wiki for more complete discussion of these options; NOTE: ALL ARGUMENTS ARE COMPULSORY AND NO DEFAULTS ARE SET
#-o = ouput dir
#-f= format of input mutations; this must be one of vcf (input mutations are as vcf files, 1 per sample, in which case they must be placed under a single directory) or bed (input mutations are in bed format: chr,start,stop,samplename, in a single file)
#-v = input mutations; if the format is vcf, then this must be a path to raw VCF files directory and all files must be placed directly under this directory; if the format is bed, then this must be the path to the bed file
#-g = genome build; must be hg19 or hg38; will determine the blacklist, annotation file, and snp file used
#-s = snp filter option; must be one of: n (none; all snps will be kept in), a (all, all snps from dbSNP 147 will be removed), c (common; snps from dbSNP 147 Common will be removed; these have a minor allele frequency >=1%); an extra option "o" is allowed when the input is in vcf format, in which the old smurf way of filtering snps is done (removing rs variants by ID)
#-r = genomic regions of interest in bed format; should have 4 columns:chr,start,end,name (for example sample name but can be anything)
#-m = method to use for BMR calculation; should be either 'global' (=bmr is total number of mutations (Across all samples) in regions divided by total bp in regions) or 'peak' (=bmr is  calculated based on which samples have mutations in that region; the individual sample bmr is calculated as for global, and the average of the BMRs for the samples represented in the peak is taken as the peak BMR.)

while getopts o:f:v:d:g:s:r:m: option
do
        case "${option}"
        in
        o) OUTDIR=${OPTARG};;
		f) FORMAT=${OPTARG};;
        v) VARIANTS=${OPTARG};;
        g) GENOME=${OPTARG};;
		s) SNP=${OPTARG};;
        r) REGIONS=${OPTARG};;
		m) METHOD=${OPTARG};;
        esac
done

module load R/3.3.0
module load bedtools/2.23.0

mkdir -p $OUTDIR


if [ $FORMAT = "vcf" ]
then
	echo Your have chosen to input vcf files
	num_samples=$(ls $VARIANTS/*vcf | wc -l)

	# Create list of files
	ls -d1 $VARIANTS/*vcf > $OUTDIR/raw_files_full_path.txt
	# remove path from file name
	touch $OUTDIR/file_list_temp.txt
	cat $OUTDIR/raw_files_full_path.txt | while read i; do
		basename $i >> $OUTDIR/file_list_temp.txt
	done
	mv $OUTDIR/file_list_temp.txt $OUTDIR/file_list.txt


	# For each VCF file, filter out variants in blacklisted regions and remove snps
	if [ $GENOME = "hg19" ]
	then
		echo The chosen genome build is hg19
		if [ $SNP = "a" ]
		then
			echo Filters: All SNPs + blacklist will be filtered out
			cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -sorted  -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg19_snp147All.bed  > $OUTDIR/$i
			done
		elif [ $SNP = "n" ]
		then
			echo Filters: blacklist will be filtered out
			cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -sorted -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed  > $OUTDIR/$i
			done
		elif [ $SNP = "c" ]
		then
			echo Filters: Common SNPs + blacklist will be filtered out
			cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -sorted -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg19_snp147Common.bed > $OUTDIR/$i
			done
		elif [ $SNP = "o" ]
                then
                        echo Filters: Legacy SNP removal by ID will be done  + blacklist will be filtered out
                        cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed > $OUTDIR/$i
				grep -F '#' $OUTDIR/$i > $OUTDIR/vcfheader.tmp
				grep -Fv '#' $OUTDIR/$i | sed -n '/rs/!p' > $OUTDIR/vcfbody.tmp
				cat $OUTDIR/vcfheader.tmp $OUTDIR/vcfbody.tmp > $OUTDIR/$i
                        done
		else
			echo Incorrect SNP argument
			exit 1
		fi
	elif [ $GENOME = "hg38" ]
	then
		echo The chosen genome build is hg38
		if [ $SNP = "a" ]
		then
			echo Filters: All SNPs + blacklist will be filtered out
			cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -sorted -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg38_snp147All.bed  > $OUTDIR/$i
			done
		elif [ $SNP = "n" ]
		then
			echo Filters: blacklist will be filtered out
			cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -sorted -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed  > $OUTDIR/$i
			done
		elif [ $SNP = "c" ]
		then
			echo Filters: Common SNPs + blacklist will be filtered out
			cat $OUTDIR/file_list.txt | while read i; do
				bedtools intersect -v -header -sorted -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg38_snp147Common.bed > $OUTDIR/$i
			done
                elif [ $SNP = "o" ]
                then
                        echo Filters: Legacy SNP removal by ID will be done  + blacklist will be filtered out
                        cat $OUTDIR/file_list.txt | while read i; do
                                bedtools intersect -v -header -a $VARIANTS/$i -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed > $OUTDIR/$i
                                grep -F '#' $OUTDIR/$i > $OUTDIR/vcfheader.tmp
                                grep -Fv '#' $OUTDIR/$i | sed -n '/rs/!p' > $OUTDIR/vcfbody.tmp
                                cat $OUTDIR/vcfheader.tmp $OUTDIR/vcfbody.tmp > $OUTDIR/$i
                        done
		else
			echo Incorrect SNP argument
			exit 1
		fi
	else
		echo Incorrect genome argument
		exit 1
	fi


 
	# Create list of VCF files with full path
	ls -d1 $OUTDIR/*vcf > $OUTDIR/filt_files_full_path.txt
	# Concatenate file paths into one string to use in bedintersect command
	myfiles=$(cat $OUTDIR/filt_files_full_path.txt | tr -s '\n' ' ')

	# Concatenante file names (without path and extension) to use in bedintersect command
	touch $OUTDIR/tmpfile.txt
	cat $OUTDIR/file_list.txt | while read i; do
	    basename $i >> $OUTDIR/tmpfile.txt
	done
	cat $OUTDIR/tmpfile.txt|sed 's/.vcf//'|sed 's/_.*//' > $OUTDIR/names.txt
	mynames=$(cat $OUTDIR/names.txt| tr -s '\n' ' ')

	#Make concatenated VCF file intersected with regions of interest
	bedtools intersect -wa -wb \
	    -a $REGIONS \
	    -b $myfiles \
	    -names $mynames | cut -f 1-7 > $OUTDIR/FiltVarsInPeaks.txt
 
	#create list of VCF sample names with number of variants in each
	touch $OUTDIR/counts_tmp.txt
	cat $OUTDIR/file_list.txt | while read i; do
	    grep -v "^#" $OUTDIR/$i | wc -l >> $OUTDIR/counts_tmp.txt
	done
	mv $OUTDIR/counts_tmp.txt $OUTDIR/counts.txt

	paste $OUTDIR/names.txt $OUTDIR/counts.txt > $OUTDIR/named_counts.txt
 
	rm $OUTDIR/tmpfile.txt
	rm $OUTDIR/{counts.txt,file_list.txt,filt_files_full_path.txt,names.txt,raw_files_full_path.txt}




elif [ $FORMAT = "bed" ]
then
	echo Your have chosen to input a bed file of variants
	#sort bed file
	sort -k1,1 -k2,2n $VARIANTS > $OUTDIR/varinput.srt.bed
	# Filter out variants in blacklisted regions and remove snps
	if [ $GENOME = "hg19" ]
	then
		echo The chosen genome build is hg19
		if [ $SNP = "a" ]
		then
			echo Filters: All SNPs + blacklist will be filtered out
			bedtools intersect -v -header -sorted -a $OUTDIR/varinput.srt.bed -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg19_snp147All.bed  > $OUTDIR/varinput.srt.filt.bed
		elif [ $SNP = "n" ]
		then
			echo Filters: blacklist will be filtered out
			bedtools intersect -v -header -sorted -a $OUTDIR/varinput.srt.bed -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed  > $OUTDIR/varinput.srt.filt.bed
		elif [ $SNP = "c" ]
		then
			echo Filters: Common SNPs + blacklist will be filtered out
			bedtools intersect -v -header -sorted -a $OUTDIR/varinput.srt.bed -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg19.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg19_snp147Common.bed > $OUTDIR/varinput.srt.filt.bed
		else
			echo Incorrect SNP argument
			exit 1
		fi
	elif [ $GENOME = "hg38" ]
	then
		echo The chosen genome build is hg38
		if [ $SNP = "a" ]
		then
			echo Filters: All SNPs + blacklist will be filtered out
			bedtools intersect -v -header -sorted -a $OUTDIR/varinput.srt.bed -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg38_snp147All.bed  > $OUTDIR/varinput.srt.filt.bed
		elif [ $SNP = "n" ]
		then
			echo Filters: blacklist will be filtered out
			bedtools intersect -v -header -sorted -a $OUTDIR/varinput.srt.bed -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed  > $OUTDIR/varinput.srt.filt.bed
		elif [ $SNP = "c" ]
		then
			echo Filters: Common SNPs + blacklist will be filtered out
			bedtools intersect -v -header -sorted -a $OUTDIR/varinput.srt.bed -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/blacklist_hg38.bed /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/hg38_snp147Common.bed > $OUTDIR/varinput.srt.filt.bed
		else
			echo Incorrect SNP argument
			exit 1
		fi
	else
		echo Incorrect genome argument
		exit 1
	fi
	
	#Make concatenated VCF file intersected with regions of interest
        bedtools intersect -wa -wb \
            -a $REGIONS \
            -b $OUTDIR/varinput.srt.filt.bed | awk -F'\t' -v OFS="\t" '{ print $1, $2, $3, $4, $8, $5, $6 }' > $OUTDIR/FiltVarsInPeaks.txt


	#create list of VCF sample names with number of variants in each
	cut -f 4 $OUTDIR/varinput.srt.filt.bed | sort|uniq -c > $OUTDIR/counts_tmp.txt #extract unique counts of each sample in variant file
	sed -e 's/^ *//' $OUTDIR/counts_tmp.txt > $OUTDIR/counts_tmp2.txt #remove whitespaces from start of line
	awk -F"[ ]" '{printf $2 "\t" $1 "\n"}' $OUTDIR/counts_tmp2.txt > $OUTDIR/named_counts.txt #print sample name and count on each line sep by tab
	rm $OUTDIR/counts_tmp*
else
	echo Incorrect format argument
	exit 1	
fi


#Annotate genomic regions with promoter or 'DistalRE' using Gencode annotation
if [ $GENOME = "hg19" ]
then
	mapBed -c 4 -o distinct -null "DistalRE" -a $REGIONS -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/Gencodev19_Promoters_NEW.bed > $OUTDIR/annotated_regions.txt
elif [ $GENOME = "hg38" ]
then
	mapBed -c 4 -o distinct -null "DistalRE" -a $REGIONS -b /mnt/work1/users/lupiengroup/CommonScripts/SMURF/annotation/Gencodev24_Promoters_NEW.bed > $OUTDIR/annotated_regions.txt
else
	echo Incorrect genome argument
	exit 1
fi	



echo Successfully filtered and annotated input files, now starting hotspot analysis

#Run hotspot analysis
Rscript /mnt/work1/users/lupiengroup/CommonScripts/SMURF/smurf.r "$OUTDIR" "$METHOD"
