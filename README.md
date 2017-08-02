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
