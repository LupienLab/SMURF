#!/bin/bash
#$ -cwd
#$ -N smurf
#$ -j y
#$ -q lupiengroup
#$ -S /bin/bash

#-o = ouput dir
#-v = path to raw VCF files; all files must be placed directly under this folder; at the moment sample names are created by removing portions of file names after "_"
#-s = separator in the names of the vcf files; the string given here and any character that follows will be removed from the vcf file names to create the sample names; eg: "sample1_cohort1.vcf", with << -s _ >> the sample name will be "sample1"; special regex characters cannot be used alone (eg ? or .) 
#-b = encode blacklist 
#-r = genomic regions of interest in bed format; should have 4 columns:chr,start,end,name(for example sample name but can be anything)
#-a = annotation file for promoters, eg:gencode
#-m = method to use for BMR calculation; should be either 'global' (=bmr is total number of mutations (Across all samples) in regions divided by total bp in regions) or 'peak' (=bmr is  calculated based on which samples have mutations in that region; the individual sample bmr is calculated as for global, and the average of the BMRs for the samples represented in the peak is taken as the peak BMR.)

while getopts o:v:s:b:r:a:m: option
do
        case "${option}"
        in
                o) OUTDIR=${OPTARG};;
                v) VCFDIR=${OPTARG};;
				s) MYSEP=${OPTARG};;
                b) BLACKLIST=${OPTARG};;
                r) REGIONS=${OPTARG};;
				a) ANNO=${OPTARG};;
				m) METHOD=${OPTARG};;	
        esac
done

module load R/3.3.0
module load bedtools/2.23.0

mkdir -p $OUTDIR

num_samples=$(ls $VCFDIR/*vcf | wc -l)

# Create list of files
ls -d1 $VCFDIR/*vcf > $OUTDIR/raw_files_full_path.txt
# remove path from file name
touch $OUTDIR/file_list_temp.txt
cat $OUTDIR/raw_files_full_path.txt | while read i; do
	basename $i >> $OUTDIR/file_list_temp.txt
done
mv $OUTDIR/file_list_temp.txt $OUTDIR/file_list.txt


# For each VCF file, filter out variants in blacklisted regions and remove snps (any variant annotated to an rs number)
cat $OUTDIR/file_list.txt | while read i; do
	bedtools intersect -v -header \
	-a $VCFDIR/$i \
	-b $BLACKLIST \
	> $OUTDIR/$i

	grep -F '#' $OUTDIR/$i > $OUTDIR/vcfheader.tmp
	grep -Fv '#' $OUTDIR/$i | sed -n '/rs/!p' > $OUTDIR/vcfbody.tmp
	cat $OUTDIR/vcfheader.tmp $OUTDIR/vcfbody.tmp > $OUTDIR/$i
done

 
# Create list of files
ls -d1 $OUTDIR/*vcf > $OUTDIR/filt_files_full_path.txt
#concatenate files into one string to use in bedintersect command
myfiles=$(cat $OUTDIR/filt_files_full_path.txt | tr -s '\n' ' ')

#concatenante file names (without path and extension) to use in bedintersect command
touch $OUTDIR/tmpfile.txt
cat $OUTDIR/file_list.txt | while read i; do
    basename $i >> $OUTDIR/tmpfile.txt
done
 
cat $OUTDIR/tmpfile.txt|sed "s/$MYSEP.*//" > $OUTDIR/names.txt
mynames=$(cat $OUTDIR/names.txt| tr -s '\n' ' ')

#Make concatenanted VCF file
bedtools intersect -wa -wb \
    -a $REGIONS \
    -b $myfiles \
    -names $mynames \
    > $OUTDIR/fullVCF.txt
 
#create list of VCF sample names with number of variants in each
touch $OUTDIR/counts_tmp.txt
cat $OUTDIR/file_list.txt | while read i; do
    grep -v "^#" $OUTDIR/$i | wc -l >> $OUTDIR/counts_tmp.txt
done
mv $OUTDIR/counts_tmp.txt $OUTDIR/counts.txt

paste $OUTDIR/names.txt $OUTDIR/counts.txt > $OUTDIR/named_counts.txt
 
rm $OUTDIR/tmpfile.txt


#Annotate genomic regions with promoter/Enhancer category
mapBed -c 4 -o distinct -null "DistalRE" -a $REGIONS -b $ANNO > $OUTDIR/annotated_regions.txt

#Run hotspot analysis
Rscript /mnt/work1/users/lupiengroup/CommonScripts/SMURF/smurf_v2.r "$OUTDIR" "$METHOD"
