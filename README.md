# SMURF
**S**ignificantly **Mu**tated **R**egion **F**inder

SMURF identifies significantly mutated genomic regions in a set of samples.

Developed by Paul Guilhamon

# Installation

Clone this repo via `git clone`.

The SNP147 annotation files from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) are very large, so they are not included in this repository.

## Requirements

SMuRF is implemented in Bash and R; it runs on any platform with Bash (≥4.1.2), R (≥3.3.0) and BEDTools (≥2.23.0).
It requires the following R packages:

* GenomicRanges
* gtools
* gplots
* ggplot2
* data.table
* psych
* dplyr

# Usage

```shell
qsub smurf.sh \
-o OUTDIR \
-f [vcf|bam] \
-v VARPATH \
-g [hg19|hg38] \
-s [n|c|a|o] \
-r REGIONS \
-m [global|peak]
```

## Arguments

| Parameter | Description |
|-----------|-------------|
| `-o` | Output directory (absolute path) |
| `-f [vcf\|bam]` | Format chosen for input variants. **vcf**: individual VCF formatted files for each sample must be placed under a single directory. The full path to that directory is what should be supplied to the `-v` argument described below. **NOTE**: when using this option the sample names used for output files will be automatically generated from the names of the VCF files themselves: if the filename contains an underscore, then the sample name will be everything that comes before the underscore; otherwise the `.vcf` at the end of the filename will be clipped off and the rest will be used as as sample name. Do make sure your sample names do not become identical once the portion of the filename after the underscore is removed. **bed**: a single BED file containing all variants to be analysed. The format must be the following: `chr start stop name`. The full path to that bed file is what should be supplied to the `-v` argument described below. |
| `-v` | The input variants themselves; depending on the format chosen above, either the path to a directory containing individual VCF files OR the path to a single BED file. |
| `-s [n\|c\|a\|o]` | SNP filter to apply to the input variants. **n**: no SNP filter, only blacklisted regions will be removed. **c**: blacklist + Common snps from dbSNP147 will be removed; *recommended option for most uses*. **a**: blacklist + all snps from dbSNP147 will be removed. Additionally, only when using the VCF files as input, the option `o` will be valid here, and will provide snp filtering as done by the previous version of SMURF: removal of any variant annotated to an rs ID within the VCF file. |
| `-p` | Promoter annotation file; those derived from Gencodev19 (human genome build hg19) and Gencodev24 (human genome build hg38) are provided in the SMuRF directory. These can be replaced with any annotation of your choice, but the output figures by SMuRF will automatically colour-code regions by either the gene name linked to a promoter region or the mention "DistalRE" for distal regulatory element, as this tool was originally designed to look for mutation-enriched non-coding regions. |
| `-r` | Genomic regions of interest that will be checked for significant enrichment of mutations. Must be supplied in the BED format specified by `-f`. The 4th column is required, but can be any string. Having the name there could make downstream analyses easier but isn't required by SMURF. |
| `-m [allsamples\|regionsamples]` | Method used to calculate the background mutation rate (BMR) in the binomial test. **allsamples**: BMR is total number of mutations (across all samples) in regions divided by total bp in regions. **regionsamples**: BMR is calculated based on which samples have mutations in that region. The individual sample BMR is calculated as for global, and the average of the BMRs for the samples represented in the region is taken as the region BMR. |

## Input

* Variants file(s): These should contain variants of interest. SMuRF has options for applying SNP filtering but will not perform arbitrary quality filters. They should be in either VCF format as multiple files (1 per sample) under the same directory OR in BED format, with all mutations for all samples within the same file, and exactly 4 columns: chr/start/stop/name.
  **NOTE**: if using the BED format, keep one row per mutation per sample, ie: if a mutation is recurrent and found in multiple samples these should be placed on separate rows. Also note that the start and stop positions are NOT the same; the start column corresponds to the POS column of a VCF-formatted file minus 1; (BED and VCF follow different numbering conventions; see https://www.biostars.org/p/84686/ for more information).
  See format example below:

  | chr | start | end | name |
  |-----|-------|-----|------|
  | chr10 | 101776719 | 101776719 | Sample1 |
  | chr10 | 101877882 | 101877882 | Sample2 |
  | ... | ... | ... | ... |

* Genomic Regions: a BED file of genomic regions of interest. The format is chr/start/stop/name ; Note: the 'name' column is required but can be any string; it is recommended to use it for the names of biological samples from which each genomic region was identified for easy trace-back in downstream analyses, although this is not required for SMuRF. These genomic regions and the variants specified above do not need to come from matched samples. Unlike the set of input variants above, it is important that the genomic regions be mutually exclusive. If in your cohort regions from two different samples overlap, they should be merged as one and the sample name column can have all the contributing sample names in a comma-separated list.
  See format example below:

| chr | start | end | name |
|-----|-------|-----|------|
| chr1 | 10027 | 10612 | Sample1 |
| chr1 | 237701 | 237901 | Sample3,Sample5,Sample22 |
| ... | ... | ... | ... |

* Annotation files: GENCODE promoter annotation, and dbSNP sets for both hg19 and hg38 are provided within the SMuRF directory.
  * Promoter annotation: The Gencode set of TSS coordinates for each genome build (hg19: v19, hg38: v24) was downloaded from the UCSC table browser. Promoters were annotated from -2.5kb to +0.5kb around the TSS.
  * SNPs: dbSNP sets v147 were downloaded from the UCSC table browser for each genome build. The "Common" version of the database is made available with the tool; this can be used to remove only those that have a minor allele frequency of >= 1%, of use if one wants to avoid removing potentially interesting somatic mutations that also appear as very rare polymorphisms in the population. The larger file containing all SNPs (including those found at frequency < 1% ) can be downloaded from the dbSNP website.

## Output

SMuRF outputs a number of files, some for display and evaluation purposes, others as end results, and a final set for input into downstream tools, such as [C3D](https://github.com/LupienLabOrganization/C3D) or [GO](http://www.geneontology.org) annotation.

### Results Files

| File | Description |
|------|-------------|
| filtered vcf files | VCF files will be filtered for SNPs and blacklist regions a specified by the user |
|`FiltVarsInPeaks.txt` | Concatenated list of unique variants, overlapped with the input genomic regions, and including the names of samples they were identified in |
| `summary.txt` | Summary of the number of input regions,number mutated, number mutated significantly, and their distribution in promoter and enhancer regions |
| `annotated_regions.txt` | All input genomic regions annotated to either a promoter region (using Gencode as described above) OR to "DistalRE" (=distal regulatory element) if they fall outside promoter regions |
| `named_counts.txt` | Count of unique variants in each sample |
| `Mutated_Regions.txt` and other files with that prefix | All information on the mutated regions and stats on significance of mutation. Other files include suffixes indicating the subsequent filtering steps. `Freq3`: only includes regions with >=3 mutations from >= 3 different samples. `qval0.05`: only includes regions significantly mutated |

### Results Plots

| File | Description |
|------|-------------|
| `QQplot_Freq3.pdf` | Self-explanatory |
| `barplot_mut_inRegions_Counts_and_Perc.pdf` | Summary plots of the number of variants per sample and how many overlap with the provided genomic intervals |
| `Mutation_Rate_Plot_Freq3_qval0.05.pdf` | Scatter plot of significantly mutated regions:-log10(q-value) vs region mutation rate |
| `Mutation_Rate_Plot_Freq3_qval0.05_Promoters_VS_DistalREs.pdf` | Same as above with added promoter/DistalRE annotation |
| `Sample_Frequency_Plot_Freq3_qval0.05.pdf` and `Sample_Frequency_Plot_Freq3_qval0.05_Promoters_VS_DistalREs.pdf` | Same as the two above, but with the number of unique samples mutated in the region as the x-axis as opposed tot he mutation rate |
| `Sample_Frequency_Plot_Freq3_qval0.05_DistalREs.pdf` and `Sample_Frequency_Plot_Freq3_qval0.05_Promoters.pdf` | Same as `Sample_Frequency_Plot_Freq3_qval0.05_Promoters_VS_DistalREs.pdf` but separating out those regions annotated to promoters and DistalREs |

### Files for use in dowstream analyses

| File | Description |
|------|-------------|
| `Mutated_Regions_sigmut0.05_DistalRE_forc3d.bed` | Significantly mutated regions annotated to enhancers, to use as input for C3D to identify downstream target promoters |
| `Mutated_Regions_sigmut0.05_DistalReAndProm.bed` | All significantly mutated regions; for use in GO applications for example, or TF motif enrichment |
| `Mutated_Regions_allmut_DistalREAndProm.bed` | All mutated regions |
