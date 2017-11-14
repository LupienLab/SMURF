#!/bin/bash


#This file is part of SMuRF

#Copyright 2017 Paul Guilhamon


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



#Brief description of arguments; please refer to README for more complete discussion of these options; NOTE: ALL ARGUMENTS ARE COMPULSORY AND NO DEFAULTS ARE SET:

# -o = ouput dir
# -f= format of input mutations; this must be one of vcf (input mutations are as vcf files, 1 per sample, in which case they must be placed under a single directory) or bed (input mutations are in bed format: chr,start,stop,samplename, in a single file)
# -v = input mutations; if the format is vcf, then this must be a path to raw VCF files directory and all files must be placed directly under this directory; the sample names used in otuput figures will be derived from the file names by simply removing the ".vcf" extension; if the format is bed, then this must be the path to the bed file
# -s = snp filter option; must be one of: n (none; all snps will be kept in), full path to snp file of choice in bed format
# -p = promoter annotation file; file containing the annotated promoters for the genome of interest
# -r = genomic regions of interest in bed format; should have 4 columns:chr,start,end,name (for example sample name but can be anything)
# -m = method to use for BMR calculation; should be either 'allsamples' (=bmr is total number of mutations (Across all samples) in regions divided by total bp in regions) or 'regionsamples' (=bmr is  calculated based on which samples have mutations in that region; the individual sample bmr is calculated as for allsamples, and the average of the BMRs for the samples represented in the region is taken as the region BMR.)




DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while getopts o:f:v:p:s:r:m: option
do
    case "${option}"
        in
        o) OUTDIR=${OPTARG};;
        f) FORMAT=${OPTARG};;
        v) VARIANTS=${OPTARG};;
        s) SNP=${OPTARG};;
        p) PROM=${OPTARG};;
        r) REGIONS=${OPTARG};;
        m) METHOD=${OPTARG};;
    esac
done

mkdir -p "$OUTDIR"


if [ "$FORMAT" = "vcf" ]
then
    echo Your have chosen to input vcf files
    num_samples=$(ls $VARIANTS/*vcf | wc -l)

    # Create list of files
    ls -d1 $VARIANTS/*vcf > "$OUTDIR/raw_files_full_path.txt"
    # remove path from file name
    touch "$OUTDIR/file_list_temp.txt"
    cat "$OUTDIR/raw_files_full_path.txt" | while read i;
    do
        basename $i >> "$OUTDIR/file_list_temp.txt"
    done
    mv "$OUTDIR/file_list_temp.txt" "$OUTDIR/file_list.txt"

    # For each VCF file,remove snps
    if [ "$SNP" = "n" ]
    then
        echo "Filters: No SNP filter applied; VCFs will be sorted only."
        cat "$OUTDIR/file_list.txt" | while read i;
        do
            grep -F '#' "$VARIANTS/$i" > "$OUTDIR/vcfheader.tmp"
            grep -Fv '#' "$VARIANTS/$i" > "$OUTDIR/vcfbody.tmp"
            sort -k1,1 -k2,2n "$OUTDIR/vcfbody.tmp" > "$OUTDIR/vcfbody2.tmp"
            cat "$OUTDIR/vcfheader.tmp" "$OUTDIR/vcfbody2.tmp" > "$OUTDIR/$i"
        done
    else
        echo "Filters: VCFs will be sorted and filtered with provided file"
        cat "$OUTDIR/file_list.txt" | while read i;
        do
            grep -F '#' "$VARIANTS/$i" > "$OUTDIR/vcfheader.tmp"
            grep -Fv '#' "$VARIANTS/$i" > "$OUTDIR/vcfbody.tmp"
            sort -k1,1 -k2,2n "$OUTDIR/vcfbody.tmp" > "$OUTDIR/vcfbody2.tmp"
            cat "$OUTDIR/vcfheader.tmp" "$OUTDIR/vcfbody2.tmp" > "$OUTDIR/tmp.vcf"
            bedtools intersect -v -header -sorted -a "$OUTDIR/tmp.vcf" -b "$SNP" > "$OUTDIR/$i"
        done
    fi

    rm "$OUTDIR/tmp.vcf"
    # Create list of VCF files with full path
    ls -d1 $OUTDIR/*vcf > "$OUTDIR/filt_files_full_path.txt"
    # Concatenate file paths into one string to use in bedintersect command
    myfiles=$(cat $OUTDIR/filt_files_full_path.txt | tr -s '\n' ' ')

    # Concatenante file names (without path and extension) to use in bedintersect command
    touch "$OUTDIR/tmpfile.txt"
    cat "$OUTDIR/file_list.txt" | while read i;
    do
        basename "$i" >> "$OUTDIR/tmpfile.txt"
    done
    cat "$OUTDIR/tmpfile.txt" | sed 's/.vcf//' > "$OUTDIR/names.txt"
    mynames=$(cat $OUTDIR/names.txt| tr -s '\n' ' ')

    #Make concatenated VCF file intersected with regions of interest
    bedtools intersect -wa -wb \
        -a "$REGIONS" \
        -b "$myfiles" \
        -names "$mynames" | cut -f 1-7 > "$OUTDIR/FiltVarsInPeaks.txt"
 
    #create list of VCF sample names with number of variants in each
    touch "$OUTDIR/counts_tmp.txt"
    cat "$OUTDIR/file_list.txt" | while read i;
    do
        grep -v "^#" "$OUTDIR/$i" | wc -l >> "$OUTDIR/counts_tmp.txt"
    done
    mv "$OUTDIR/counts_tmp.txt" "$OUTDIR/counts.txt"

    paste "$OUTDIR/names.txt" "$OUTDIR/counts.txt" > "$OUTDIR/named_counts.txt"
 
    rm "$OUTDIR/tmpfile.txt"
    rm $OUTDIR/{counts.txt,file_list.txt,filt_files_full_path.txt,names.txt,raw_files_full_path.txt}

elif [ "$FORMAT" = "bed" ]
then
    echo "You have chosen to input a bed file of variants"
    #sort bed file
    sort -k1,1 -k2,2n "$VARIANTS" > "$OUTDIR/varinput.srt.bed"
    # Filter out snps
    if [ "$SNP" = "n" ]
    then
        echo "Filters: No SNP filter applied"
        cp "$OUTDIR/varinput.srt.bed" "$OUTDIR/varinput.srt.filt.bed"
    else
        echo "Filters: Variants will be filtered with provided file"
        bedtools intersect -v -header -sorted \
            -a "$OUTDIR/varinput.srt.bed" \
            -b "$SNP" \
            > "$OUTDIR/varinput.srt.filt.bed"
    fi

    #Make concatenated VCF file intersected with regions of interest
    bedtools intersect -wa -wb \
        -a "$REGIONS" \
        -b "$OUTDIR/varinput.srt.filt.bed" \
        | awk -F '\t' -v OFS="\t" '{ print $1, $2, $3, $4, $8, $5, $6 }' \
        > "$OUTDIR/FiltVarsInPeaks.txt"

    #create list of VCF sample names with number of variants in each
    cut -f 4 "$OUTDIR/varinput.srt.filt.bed" | sort | uniq -c > "$OUTDIR/counts_tmp.txt"
    #extract unique counts of each sample in variant file
    sed -e 's/^ *//' "$OUTDIR/counts_tmp.txt" > "$OUTDIR/counts_tmp2.txt"
    #remove whitespaces from start of line
    awk -F "[ ]" '{printf $2 "\t" $1 "\n"}' "$OUTDIR/counts_tmp2.txt" > "$OUTDIR/named_counts.txt"
    #print sample name and count on each line sep by tab
    rm $OUTDIR/counts_tmp*
else
    echo "Incorrect format argument"
    exit 1
fi

#Annotate genomic regions with promoter or 'DistalRE' using Gencode annotation
mapBed -c 4 -o distinct -null "DistalRE" -a "$REGIONS" -b "$PROM" > "$OUTDIR/annotated_regions.txt"

echo "Successfully filtered and annotated input files, now starting hotspot analysis"

#Run hotspot analysis
Rscript "$DIR/smurf.r" "$OUTDIR" "$METHOD"