# SMURF
**S**ignificantly **Mu**tated **R**egion **F**inder

SMURF identifies significantly mutated genomic regions in a set of samples

Developed by Paul Guilhamon

# Installation

Clone this repo via `git clone`.

The SNP147 annotation files from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) are very large, so they are not included in this repository.

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

## Parameters

| Parameter | Description |
|-----------|-------------|
| `-o` | Output directory (absolute path) |
| `-f [vcf\|bam]` | Format chosen for input variants. **vcf**: individual VCF formatted files for each sample must be placed under a single directory. The full path to that directory is what should be supplied to the `-v` argument described below. **NOTE**: when using this option the sample names used for output files will be automatically generated from the names of the VCF files themselves: if the filename contains an underscore, then the sample name will be everything that comes before the underscore; otherwise the `.vcf` at the end of the filename will be clipped off and the rest will be used as as sample name. Do make sure your sample names do not become identical once the portion of the filename after the underscore is removed. **bed**: a single BED file containing all variants to be analysed. The format must be the following: `chr start stop name`. The full path to that bed file is what should be supplied to the `-v` argument described below. |
| `-v` | The input variants themselves; depending on the format chosen above, either the path to a directory containing individual VCF files OR the path to a single BED file. |
| `-g [hg19\|hg38]` | Genome build. This will determine the build of annotation files used by SMURF. Obviously, the input variants and input genomic regions should both be of the chosen build. |
| `-s [n\|c\|a\|o]` | SNP filter to apply to the input variants. **n**: no SNP filter, only blacklisted regions will be removed. **c**: blacklist + Common snps from dbSNP147 will be removed; *recommended option for most uses*. **a**: blacklist + all snps from dbSNP147 will be removed. Additionally, only when using the VCF files as input, the option `o` will be valid here, and will provide snp filtering as done by the previous version of SMURF: removal of any variant annotated to an rs ID within the VCF file. |
| `-r` | Genomic regions of interest that will be checked for significant enrichment of mutations. Must be supplied in the BED format specified by `-f`. The 4th column is required, but can be any string. Having the name there could make downstream analyses easier but isn't required by SMURF. |
| `-m [global\|peak]` | Method used to calculate the background mutation rate (BMR) in the binomial test. **global**: BMR is total number of mutations (across all samples) in regions divided by total bp in regions. **peak**: BMR is calculated based on which samples have mutations in that region. The individual sample BMR is calculated as for global, and the average of the BMRs for the samples represented in the peak is taken as the peak BMR. |

## Input

* Variants file(s): These should contain variants of interest. SMURF has options for applying SNP and blacklist filtering but will not perform arbitrary quality filters. They should be in either VCF format as multiple files (1 per sample) under the same directory OR in BED format, with all mutations for all samples within the same file, and exactly 4 columns: chr/start/stop/name. **NOTE**: if using the BED format, keep one row per mutation per sample, ie: if a mutation is recurrent and found in multiple samples these should be placed on separate rows. Also note that the start and stop positions are the same; these correspond to the POS column of a VCF-formatted file (BED and VCF follow different numbering conventions). See format example below:

| chr | start | end | name |
|-----|-------|-----|------|
| chr10 | 101776719 | 101776719 | Sample1 |
| chr10 | 101877882 | 101877882 | Sample2 |
| ... | ... | ... | ... |

* Genomic Regions: a BED file of genomic regions of interest. The format is `chr start stop name`. **Note**: the `name` column is required but can be any string; it is recommended to use it for the names of biological samples from which each genomic region was identified for easy traceback in downstream analyses, although this is not required for SMURF. These genomic regions and the variants specified above do not need to come from matched samples. Unlike the set of input variants above, it is important that the genomic regions be mutually exclusive. If in your cohort regions from two different samples overlap, they should be merged as one and the sample name column can have all the contributing sample names in a comma-separated list. See format example below:

| chr | start | end | name |
|-----|-------|-----|------|
| chr1 | 10027 | 10612 | Sample1 |
| chr1 | 237701 | 237901 | Sample3,Sample5,Sample22 |
| ... | ... | ... | ... |

* Annotation files: ENCODE blacklist, Gencode promoter annotation, and dbSNP sets (All and Common) for both hg19 and hg38 are provided within the annotation directory. The user on the cluster does not need to provide any of these.
  * [Blacklist](https://sites.google.com/site/anshulkundaje/projects/blacklists) generated by Anshul Kundaje. They represent "artefact regions that tend to show artificially high signal (excessive unstructured anomalous reads mapping)".
  * Promoter annotation: The Gencode set of TSS coordinates for each genome build (hg19: v19, hg38: v24) was downloaded from the UCSC table browser. Promoters were annotated from -2.5kb to +0.5kb around the TSS.
  * SNPs: dbSNP sets v147 were downloaded from the UCSC table browser for each genome build. Two versions of the databases are available, "All" to filter out all known SNPs and "Common" to remove only those that have a minor allele frequency of >= 1%. This latter set might be of use if one wants to avoid removing potentially interesting somatic mutations that also appear as very rare polymorphisms in the population.

## Output

SMURF outputs a number of files, some for display and evaluation purposes, others as end results, and a final set for input into downstream tools, such as [C3D](https://github.com/LupienLabOrganization/C3D) or [GO](http://www.geneontology.org) annotation.

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
| `QQplot_Freq3.pdf` | self-explanatory |
| `barplot_mut_inRegions_Counts_and_Perc.pdf` | summary plots of the number of variants per sample and how many overlap with the provided genomic intervals |
| `Mutated_Regions_Freq3_qval0.05.pdf` | scatter plot of significantly mutated regions |
| `Mutated_Regions_Freq3_qval0.05_Promoters_VS_DistalREs.pdf` | same as above with added promoter/DistalRE annotation |

### Files for use in dowstream analyses

| File | Description |
|------|-------------|
| `Mutated_Regions_sigmut0.05_DistalRE_forc3d.bed` | significantly mutated regions annotated to enhancers, to use as input for C3D to identify downstream target promoters |
| `Mutated_Regions_sigmut0.05_DistalReAndProm.bed` | all significantly mutated regions; for use in GO applications for example, or TF motif enrichment |
| `Mutated_Regions_allmut_DistalREAndProm.bed` | all mutated regions |
