#!/bin/bash
#$ -cwd
#$ -N smurf
#$ -o logs/smurf.log
#$ -j y
#$ -q lupiengroup
#$ -S /bin/bash

#-o = ouput dir
#-v = path to raw VCF files; all files must be placed directly under this folder; at the moment sample names are created by removing portions of file names after "_"
#-b = encode blacklist 
#-r = genomic regions of interest in bed format; should have 4 columns:chr,start,end,name(for example sample name but can be anything)
#-a = annotation file for promoters, eg:gencode

while getopts o:v:b:r:a: option
do
        case "${option}"
        in
                o) OUTDIR=${OPTARG};;
                v) VCFDIR=${OPTARG};;
                b) BLACKLIST=${OPTARG};;
                r) REGIONS=${OPTARG};;
		a) ANNO=${OPTARG};;
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
 
cat $OUTDIR/tmpfile.txt|sed 's/_.*//' > $OUTDIR/names.txt
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
mapBed -c 4 -o distinct -null "Enhancer" -a $REGIONS -b $ANNO > $OUTDIR/annotated_regions.txt

#Run hotspot analysis
Rscript /mnt/work1/users/lupiengroup/CommonScripts/SMURF/smurf.r "$OUTDIR"
