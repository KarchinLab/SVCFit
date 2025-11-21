
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVCFit

<!-- badges: start -->
<!-- badges: end -->

SVCFit is a fast and scalable computational tool developed to estimate
the structural variant cellular fraction (SVCF) of inversions, deletions
and tandem duplications. SVCFit is designed to run in an R environment.

## Installation

To install SVCFit from GitHub, you must have a GitHub acount. If you
don’t have an account, first sign up at [GitHub](https://github.com/).
Then, You can then install SVCFit within an R environment.

``` r
install.packages("usethis")
install.packages("remotes")

usethis::use_git_config(user.name = "Github_user_name")
usethis::create_github_token()

# This will open a web page in your browser where you can sign in to GitHub.
# Once signed in, you will be directed to a page to generate a new Personal Access Token (PAT).
# Enter a descriptive note for your PAT. Use the default "Scope Options" if you're unsure about them
# Then, Click the "Generate Token" button, and copy the generated PAT. 
# Make sure to store it securely in a text file or a password manager.

credentials::set_github_pat()

#If this is your first PAT, a pop-up screen will appear with a box to enter your PAT.
#Go ahead and enter it. Otherwise, proceed to the next step.

remotes::install_github("KarchinLab/SVCFit", build_vignettes = TRUE, dependencies = TRUE)
```

## Input your structural variants into SVCFit

SVCFit is designed to take input from Variant Call Format (VCF) files.
By default, it accepts the VCF format produced by the Manta package
\[1\]. If you have a VCF output from a structural variant caller other
than Manta, please modify to match this format:

| CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | normal | tumor |
|----|----|----|----|----|----|----|----|----|----|----|
| chr1 | 1000 | INV:6:0:1:0:0:0 | T | <INV> | . | PASS | END=1500;SVTYPE=INV;SVLEN=500 | PR:SR | 20,30:19,27 | 23,0:17,0 |
| chr2 | 5000 | DEL:7:0:1:0:0:0 | G | <DEL> | . | PASS | END=5300;SVTYPE=DEL;SVLEN=300 | PR | 15,30 | 19,0 |

## General workflow

### 1. Extract SV and SNP Information — `extract_info()`

This step loads and preprocesses the input VCF and CNV files. `extract_info()` internally performs:

1.  **Load input data** — `load_data()`
2.  **Process BND events** — `proc_bnd()`
3.  **Parse SV metadata** — `parse_sv_info()`
4.  **Parse heterozygous SNPs** — `parse_het_snp()`, `parse_snp_on_sv()`

``` r
# all in one
info <- extract_info(p_het, p_onsv, p_sv, p_cnv, flank_del=50, QUAL_tresh=100, min_alt=2, tumor_only=FALSE)

# step by step
data      <- load_data(samp, exper, chr=chr_lst, tumor_only=FALSE)
bnd_info   <- proc_bnd(data$sv, flank_del=50)
sv_info   <- parse_sv_info(data$sv, bnd_info$bnd, bnd_info$del, QUAL_tresh=100, min_alt=2)
snp_df    <- parse_het_snps(data$het_snp)
sv_phase  <- parse_snp_on_sv(data$het_on_sv, snp_df)
```

#### Function Arguments

#### `load_data()`

| Argument | Type | Default | Description |
|-----------------|-----------------|-----------------|---------------------|
| `p_het` | Character | — | Path to VCF of heterozygous SNPs. |
| `p_onsv` | Character | — | Path to VCF of SNPs overlapping SV-supporting reads. |
| `p_sv` | Character | — | Path to SV VCF. |
| `p_cnv` | Character | — | Path to CNV file. |
| `chr` | Character | — | Chromosomes to include. |
| `tumor_only` | Logical | FALSE | Whether SVs come from tumor-only calling. |

#### `proc_bnd()`

| Argument | Type | Default | Description |
|-----------------|-----------------|-----------------|---------------------|
| `sv` | data.frame | — | SV VCF data. |
| `flank_del` | numeric | 50 | Max distance to consider deletion overlapping a BND. |

#### `parse_sv_info()`

| Argument     | Type       | Default | Description                 |
|--------------|------------|---------|-----------------------------|
| `sv`         | data.frame | —       | SV VCF data.                |
| `bnd`        | data.frame | —       | Translocations (BNDs).      |
| `del`        | data.frame | —       | Deletions overlapping BNDs. |
| `QUAL_tresh` | numeric    | —       | Minimum QUAL score.         |
| `min_alt`    | numeric    | —       | Minimum alternative reads.  |

#### `parse_het_snp()`

| Argument  | Type       | Default | Description        |
|-----------|------------|---------|--------------------|
| `het_snp` | data.frame | —       | Heterozygous SNPs. |

#### `parse_snp_on_sv()`

| Argument    | Type       | Default | Description                  |
|-------------|------------|---------|------------------------------|
| `het_on_sv` | data.frame | —       | SNPs on SV-supporting reads. |
| `snp_df`    | data.frame | —       | Parsed heterozygous SNPs.    |

**Output:** A list of data frames containing parsed SV + SNP information.

### 2. Annotate SVs Using CNV and SNP Information — `characterize_sv()`

This step integrates CNV and heterozygous SNPs to infer additional information for SVs, 
including phasing, zygosity, and overlapping CNV. `characterize_sv()` internally performs:

1.  **Assign SV IDs to SNPs** — `assign_svids()`
2.  **Summarizes phasing + zygosity** — `sum_sv_info()`
3.  **Assign CNV to SV** — `assign_cnv()`
4.  **Annotate overlapping CNV** — `annotate_cnv()`, `parse_snp_on_sv()`

``` r
# all in one
sv_char <- characterize_sv(sv_phase, sv_info, cnv)

# step by step 
assign_id <- assign_svids(sv_phase, sv_info, flank=500)
sv_sum    <- sum_sv_info(sv_phase, assign_id, sv_info)
sv_cnv    <- assign_cnv(sv_sum, cnv)
anno_sv_cnv     <- annotate_cnv(sv_cnv)
```

#### Function Arguments

#### `assign_svids()`

| Argument   | Type       | Default | Description                 |
|------------|------------|---------|-----------------------------|
| `sv_phase` | data.frame | —       | Phasing/zygosity from SNPs. |
| `sv_info`  | data.frame | —       | Parsed SV metadata.         |
| `flank`    | numeric    | —       | Max assignment distance.    |

#### `sum_sv_info()`

| Argument    | Type       | Default | Description          |
|-------------|------------|---------|----------------------|
| `sv_phase`  | data.frame | —       | Phasing information. |
| `assign_id` | data.frame | —       | SNP→SV assignment.   |
| `sv_info`   | data.frame | —       | Parsed SV metadata.  |

#### `assign_cnv()`

| Argument | Type       | Default | Description     |
|----------|------------|---------|-----------------|
| `sv_sum` | data.frame | —       | Summarized SVs. |
| `cnv`    | data.frame | —       | CNV data.       |

#### `annotate_cnv()`

| Argument | Type       | Default | Description              |
|----------|------------|---------|--------------------------|
| `sv_cnv` | data.frame | —       | SVs annotated with CNVs. |

### 3. Calculate SVCF for Structural Variants — `calculate_svcf()`

This step computes the **Structural Variant Cellular Fraction (SVCF)**.

``` r
output <- calculate_svcf(input=checked, tumor_only=FALSE)
```

#### Function Arguments

| Argument      | Type       | Default | Description                            |
|--------------|--------------|--------------|--------------------------------|
| `anno_sv_cnv` | data.frame | —       | CNV-annotated SVs.                     |
| `sv_info`     | data.frame | —       | Parsed SV info.                        |
| `thresh`      | numeric    | 0.1     | Threshold for SV-before-CNV inference. |
| `samp`        | character  | —       | Sample name.                           |
| `exper`       | character  | —       | Experiment name.                       |

The output is an annotated VCF with additional fields for VAF, Rbar, r
and SVCF. VAF=variant allele frequency; Rbar=average break interval
count in a sample; r = inferred integer copy number of break intervals;
SVCF=structural variant cellular fraction.

### 4. Additional functions

*read_clone* and *attach_clone* are functions to assign structural
variants to tumor clones, when the assignment is known.

For example, if you have run a simulation and know the clonal assignment
of each structural variant, the first step is to read in your clonal
assignments.

4.1 read_clone

This function has 2 arguments:

| Argument | Type | Default | Description |
|----|----|----|----|
| `truth_path` | Character | N/A | Path to BED files storing true structural variant information with clonal assignment. Each BED file should be named like `"c1.bed, c2.bed"`, etc. Structural variants should be saved in separate BED files if they belong to different (sub)clones. |
| `mode` | Character | `inherited` | <br>- **`"inherited"`**: BED files for the child clones contain all ancestral structural variants.<br>- **`"distinct"`**: Child clones do not contain any ancestral structural variants. |

``` r
truth <- read_clone(truth_path, mode="inherited")
```
Organize your BED files according to this structure:

``` r
root/
├── true_clone/
│   ├── c1.bed/
│   ├── c2.bed/
│   ├── c3.bed/
│   └── .../
```

4.2 attach_clone

This function has 3 arguments:

| Variable | Type | Default | Description |
|----|----|----|----|
| `dat` | DataFrame | N/A | Stores structural variants for clone assignment. |
| `truth` | DataFrame | N/A | Stores the clone assignment for each structural variant. |
| `tolerance` | Integer | `6` | Sets the threshold for the maximum distance between structural variants to be considered as a single structural variant when assigning clones. |

``` r
attched <- attach_clone(dat, truth, tolerance = 6)
```
## Test datasets

Test datasets are available in this repository
https://github.com/KarchinLab/SVCFit/tree/main/inst/extdata

## Tutorial

The tutorial uses the test datasets and provides detailed instructions about how to use 
the input and the expected output.

``` r
library(SVCFit)
vignette("SVCFit_guide", package = "SVCFit")
```

## Reproducing figures in Liu et al. SVCFit: Inferring Structural Variant Cellular Fraction in Tumors.

All data used for generating the figures in Liu et al can be acquired as described below. All necessary code is available in <https://github.com/KarchinLab/SVCFit/tree/main/Paper>.

* The simulated data described in `Methods/VISOR Benchmark Set` and `Methods/Timing Comparison` is available [here](https://doi.org/10.17632/bwzb6n3xbc.1).
* The data described in `Methods/Benchmark on prostate metastasis mixtures` can be generated by accessing protected data from European Genome-phenome Archive (EGAD00001001343) and running this script
<https://github.com/mcmero/SVclone_Rmarkdown/blob/master/make_insilico_mixtures.sh> \[2\].

## Reference

1. Chen, X. et al. (2016) Manta: rapid detection of structural variants
    and indels for germline and cancer sequencing applications.
    Bioinformatics, 32, 1220-1222. <doi:10.1093/bioinformatics/btv710>
2.  Cmero, Marek, Yuan, Ke, Ong, Cheng Soon, Schröder, Jan, Corcoran,
    Niall M., Papenfuss, Tony, et al., “Inferring Structural Variant
    Cancer Cell Fraction,” Nature Communications, 11(1) (2020), 730.
