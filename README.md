
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVCFit

<!-- badges: start -->
<!-- badges: end -->

SVCFit is a super fast computational tool developed to calculate
structural variant cellular fraction (SVCF) to aid the construction of
tumor phylogeny. It takes unique characteristic of structural variants
when computing SVCF to achieve more accurate result. Currently, our
package work with Inversion, Deletion, and Tandem duplication.

## Installation

You can install the development version of SVCFit from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("KarchinLab/SVCFit")
```

## Key data structure

The input of this package should be in Variant Call Format (VCF). Each
row describes a structural variants and each column are the attributes
for this structural variants. In this example below, we have a sample
VCF from Manta (VCF from other tools need to be modified to fit this
format):

<!-- ```{r,echo = FALSE}  example_vcf=data.frame(CHROM=c("chr1","chr2"), POS=c(1000, 5000), ID=c("MantaINV:6:0:1:0:0:0","MantaDEL:7:0:1:0:0:0"), REF=c("T","G"), ALT=c("<INV>","<DEL>"), QUAL=c(".","."), FILTER=c("PASS","PASS"),INFO=c("END=1500;SVTYPE=INV;SVLEN=500","END=5300;SVTYPE=DEL;SVLEN=300"), FORMAT=c("PR:SR","PR"),tumor=c("20,30:19,27", "15,30")) # example_vcf #` -->
<!-- ``` -->

    CHROM   POS   ID    REF   ALT   QUAL    FILTER    INFO    FORMAT    tumor
    chr1    1000    MantaINV:6:0:1:0:0:0    T   <INV>   .   PASS    END=1500;SVTYPE=INV;SVLEN=500   PR:SR   20,30:19,27
    chr2    5000    MantaDEL:7:0:1:0:0:0    G   <DEL>   .   PASS    END=5300;SVTYPE=DEL;SVLEN=300   PR   15,30

## General workflow

*SVCF()* is the main function in this package that wraps all usage
described below. All functions can be ran separately.

### 1. Extract information from original VCF (extract_info)

This step assigns column names and extracts key information for
downstream calculation and filtering. Key information includes reads,
structural variant length, and structural variant coordinates.

This function has three inputs: 1. vcf_path: Character variable of path
to vcf files 2. tumor_only: Boolean variable of whether the VCF is
created without matched normal sample 3. length_filter: Numeric variable
of structural variant length filter threshold.

``` r
vcf=extract_info("~/path/to/file.vcf", tumor_only=TRUE, length_filter=0)
```

The output from *extract_info()* will be in annotated VCF format.

### 2. Check overlapping structural variants

This step checks if the genomic coordinates of two structural variants
are close enough to be considered as one event.

This function has 4 inputs: 1. dat: a dataframe to be compared 2.
compare: a dataframe used as baseline for comparison 3. tolerance: an
integer variable setting the threshold for genomic coordinates
difference to be considered as the same event 4. window: an integer
variable setting threshold for how many svructural variant to be
compared at once

Note: When *dat* and *compare* are the same dataframe, this function
will remove overlapping structural variant and use averaged reads to
represent this structural variant. When *compare* is the ground truth
from a simulation, this function will additionally remove structural
variants that are not in the simulation.

``` r
checked=check_overlap(vcf, vcf)
```

### 3. Calculate SVCF for structural variants

This step calculate structural variant cellular fraction (SVCF). The
*inferred_icn* is not an accurate estimation for integer copy number.

``` r
output <- checked %>%
      dplyr::filter(!classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = alt/(alt+0.5*ref),
             tcn= ifelse(classification =="DUP", (4*alt+2*ref)/ref, 2),
             inferred_icn = ifelse(classification=="DUP", round(tcn*vaf+2),2),
             svcf = ifelse(inferred_icn<4, tcn*vaf, tcn*vaf/(inferred_icn-2)))%>%
      dplyr::select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,tumor,classification,pos2,vaf,tcn,inferred_icn,svcf)
```

The output is an annotated VCF with additional information supporting
the calculation of SVCF.

### 4. Additional functions

*attach_clone* and *read_clone* are functions to assign SV to clones
when ground truth from simulations is known.

## Tutorial

``` r
library(SVCFit)
vignette("SVCFit_guide", package = "SVCFit")
```
