
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVCFit

<!-- badges: start -->
<!-- badges: end -->

SVCFit is a fast computational tool developed to estimate the 
structural variant cellular fraction (SVCF) of inversions, deletions and
tandem duplications. The SVCF can be used to incorporate these structural 
variants into a tumor evolutionary tree. 

## Installation

You can install SVCFit from
[GitHub](https://github.com/) with:

``` r
#set config
usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")

#Command below will generate link to github page to generate token 
usethis::create_github_token() 

#paste your PAT into pop-up that follows...
credentials::set_github_pat()

#now remotes::install_github() will work
remotes::install_github("KarchinLab/SVCFit")
```

## Input your structural variants into SVCFit

The input should be in Variant Call Format (VCF) as output by the Manta package (cite). 
If you have a VCF output from a structural variant caller other than Manta, you can modify it to match Manta format.

<!-- ```{r,echo = FALSE}  example_vcf=data.frame(CHROM=c("chr1","chr2"), POS=c(1000, 5000), ID=c("MantaINV:6:0:1:0:0:0","MantaDEL:7:0:1:0:0:0"), REF=c("T","G"), ALT=c("<INV>","<DEL>"), QUAL=c(".","."), FILTER=c("PASS","PASS"),INFO=c("END=1500;SVTYPE=INV;SVLEN=500","END=5300;SVTYPE=DEL;SVLEN=300"), FORMAT=c("PR:SR","PR"),tumor=c("20,30:19,27", "15,30")) # example_vcf #` -->
<!-- ``` -->

    CHROM   POS   ID    REF   ALT   QUAL    FILTER    INFO    FORMAT    tumor
    chr1    1000    MantaINV:6:0:1:0:0:0    T   <INV>   .   PASS    END=1500;SVTYPE=INV;SVLEN=500   PR:SR   20,30:19,27
    chr2    5000    MantaDEL:7:0:1:0:0:0    G   <DEL>   .   PASS    END=5300;SVTYPE=DEL;SVLEN=300   PR   15,30

## General workflow

*SVCF()* is the main function in this package that wraps all functionality 
described below. All functions can also be run separately. The steps executed by
*SVCF()* are:

### 1. Extract information from input VCF (extract_info)

This step assigns column names and extracts key information for
downstream calculation and filtering. Key information includes reads,
structural variant length, and structural variant coordinates.

This function has three inputs:

1.  vcf_path: Character variable of path to vcf files

2.  tumor_only: Boolean variable of whether the VCF is created without
    matched normal sample

3.  length_threshold: Numeric variable of the structural variant length filter
    threshold.

For example, the following command will generate an annotated VCF file with all structural variants with length>0
``` r
vcf=extract_info("~/path/to/file.vcf", tumor_only=TRUE, length_threshold=0)
```

The output from *extract_info()* will be in annotated VCF format.

### 2. Check overlapping structural variants

This step checks if structural variants are close enough to be considered as a single structural variant.

This function has 4 inputs:

1.  dat: a dataframe to be compared

2.  compare: a dataframe used as reference for comparison

3.  tolerance: an integer variable setting the threshold for the maximum 
    distance between structural variants to be considered as a single structural variant.

5.  window: an integer variable setting the threshold for how many
    structural variants should be checked for whether they overlap to form a single structural variant.

Note: When *dat* and *compare* are the same dataframe, this function
will merge overlapping structural variants and recompute the reads supporting the
new structural variant by taking the average of the overlapping 
structural variants. When *compare* is the ground truth
from a simulation, this function will also remove false positive structural
variants that were not included in the simulation.

``` r
checked=check_overlap(vcf, vcf)
```

### 3. Calculate SVCF for structural variants

This step calculates the structural variant cellular fraction (SVCF) for all
structural variants in the input VCF file. 

``` r
output <- checked %>%
      dplyr::filter(!classification%in%c("INS","BND"))%>%
      dplyr::mutate(vaf = alt/(alt+0.5*ref),
             tcn= ifelse(classification =="DUP", (4*alt+2*ref)/ref, 2),
             inferred_icn = ifelse(classification=="DUP", round(tcn*vaf+2),2),
             svcf = ifelse(inferred_icn<4, tcn*vaf, tcn*vaf/(inferred_icn-2)))%>%
      dplyr::select(sample, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,tumor,classification,pos2,vaf,tcn,inferred_icn,svcf)
```

The output is an annotated VCF with additional fields for VAF, Rbar, inferred ICN and SVCF.
VAF=variant allele frequency; Rbar=average break interval count in a sample; inferred ICN = inferred integer copy number; SVCF=structural variant cellular fraction

### 4. Additional functions

*attach_clone* and *read_clone* are functions to assign structural variants to tumor clones, when the assignment is known.

## Tutorial

``` r
library(SVCFit)
vignette("SVCFit_guide", package = "SVCFit")
```
