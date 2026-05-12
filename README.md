
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SVCFit <img src="inst/extdata/svcfit_logo.png" align="left" width="60"/>

SVCFit is a fast and scalable computational tool designed to estimate
the Structural Variant Cellular Fraction (SVCF) of inversions,
deletions, tandem duplications, and translocations. Developed for the R
environment, SVCFit integrates structural variant (SV) calls with Copy
Number Variation (CNV) and Single Nucleotide Polymorphism (SNP) data to
provide accurate cellular fraction estimates.

**Resources**

- Open access data: It is available on mendeley (doi:
  10.17632/2nhhdjx225.3)

- Protected Data: Available via European Genome-phenome Archive
  (EGAD00001001343).

- Prostate mixture scripts: [GitHub
  Repository](https://github.com/mcmero/SVclone_Rmarkdown/blob/master/make_insilico_mixtures.sh)

## Installation

SVCFit is hosted on GitHub. You can install it directly within R using
the `remotes` package.

**Note:** Installation requires a GitHub Personal Access Token (PAT)
because the repository is hosted on GitHub.\\

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# 1. Setup GitHub Credentials (if not already configured)
if (!requireNamespace("usethis", quietly = TRUE))
    install.packages("usethis")

# Create a token in your browser
usethis::create_github_token() 

# Store the token (paste when prompted)
credentials::set_github_pat()

# 2. Install SVCFit
remotes::install_github("KarchinLab/SVCFit", build_vignettes = TRUE, dependencies = TRUE)
```

## Input Requirements

### 1. Structural variants

SVCFit accepts standard Variant Call Format (VCF) files. By default, the
parser is optimized for VCFs produced by the SVTyper package \[2\].

| CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | normal | tumor |
|----|----|----|----|----|----|----|----|----|----|----|
| chr1 | 1000 | INV:6:0:1:0:0:0 | T | <INV> | 100 | PASS | END=1500;SVTYPE=INV;SVLEN=500;… | GT:PR:SR:… | 0/1:76,0:70,0:… | 0/1:76,0:70,0:… |
| chr2 | 5000 | DEL:7:0:1:0:0:0 | G | <DEL> | 100 | PASS | END=5300;SVTYPE=DEL;SVLEN=300;… | GT:PR:SR:… | 0/1:76,0:70,0:… | 0/1:76,0:70,0:… |

Required **INFO** fields:

- `SVTYPE` (e.g., INV, DEL, DUP, BND)
- `END`

### 2. Copy number variants

SVCFit currently utilizes copy number calls from the FACETS package
\[3\]. The tool uses allele-specific copy number, total copy number, and
the cellular fraction (cncf) to annotate SVs with overlapping CNVs. No
modification is required for standard FACETS output.

### 3. Heterozygous SNP near SV

SVCFit requires heterozygous SNP calls to phase overlapping CNVs. These
can be generated using GATK4 \[4\] HaplotypeCaller and filtered using
bcftools \[5\].

For computational efficiency, VCF can be filtered to include only
heterozygous SNPs within 500bp of SV breakpoints.

The following code is not included in SVCFit and should be run
separately.

    # 1. Call SNPs using GATK
    gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R $ref \
        -I $normal_bam \
        -O $snp_dir/SNP.vcf.gz

    # 2. Filter for heterozygous SNPs using bcftools
    bcftools view -v snps -g het -Oz -o $snp_dir/het_snp.vcf.gz $snp_dir/SNP.vcf.gz
    tabix -p vcf $snp_dir/het_snp.vcf.gz

### 4. Heterozygous SNP on SV supporting reads

To infer SV phasing, SVCFit specifically examines heterozygous SNPs
found on reads that support the structural variant. This is done
following the steps: 1) Extract SV-supporting reads from the BAM file
using samtools \[6\]. 2) Generate a read pileup using bcftools
restricted to the heterozygous SNP positions identified in step 3.

The following code is not included in SVCFit and should be run
separately.

    # Extract SV supporting reads
    samtools view -f 1 -F 2 -b $tumor_bam > $snp_dir/sup_$samp_name.bam
    samtools index $snp_dir/sup_$samp_name.bam

    # Generate pileup at known heterozygous sites
    # Note: pos_$samp_name.bed should contain the positions from het_snp.vcf.gz above
    bcftools mpileup -f $ref -a DP,AD -A \
        -R $snp_dir/pos_$samp_name.bed \
        $snp_dir/sup_$samp_name.bam -Ov > $snp_dir/het_on_sv_$samp_name.vcf

## Example Preprocessing Pipeline

The commands below provide a complete worked example for generating SVCFit inputs from raw BAM files. Each tumor-versus-normal biopsy pair is processed as a separate sample. Tool versions used in the manuscript are listed in the [Tool Versions](#tool-versions) table.

The expected input layout:

    sample/
    ├── tumor.bam        (+ .bai)
    ├── normal.bam       (+ .bai)
    ├── ref.fa           (+ .fai, .dict)
    └── intervals.bed    (optional)

### 1. Trim and align reads (if starting from FASTQ)

    trim_galore --paired tumor_R1.fq.gz tumor_R2.fq.gz -o trimmed/
    bwa mem -t 8 ref.fa trimmed/tumor_R1_val_1.fq.gz trimmed/tumor_R2_val_2.fq.gz \
      | samtools sort -@ 4 -o tumor.bam -
    samtools index tumor.bam
    # repeat for the matched normal

### 2. Mark duplicates and recalibrate (GATK4)

    gatk MarkDuplicates -I tumor.bam -O tumor.md.bam -M tumor.metrics.txt
    gatk BaseRecalibrator -I tumor.md.bam -R ref.fa --known-sites known.vcf.gz -O tumor.bqsr.table
    gatk ApplyBQSR        -I tumor.md.bam -R ref.fa --bqsr-recal-file tumor.bqsr.table -O tumor.recal.bam
    # repeat for the matched normal

### 3. Somatic SV calling (Manta)

    configManta.py --tumorBam tumor.recal.bam --normalBam normal.recal.bam \
                   --referenceFasta ref.fa --runDir manta_run/
    manta_run/runWorkflow.py
    # output: manta_run/results/variants/somaticSV.vcf.gz

For multi-caller consensus (recommended for clinical samples), run Delly and GRIDSS in parallel and merge with SURVIVOR:

    SURVIVOR merge vcf_list.txt 500 1 1 1 0 30 consensus.vcf

### 4. SV genotyping (SVtyper)

    # Ensure CIPOS and CIEND INFO fields are present; SVtyper requires them.
    # For VCFs from GRIDSS / Delly / merged consensus, run:
    # python add_cipos_ciend.py consensus.vcf > consensus.cici.vcf

    svtyper -B tumor.recal.bam -i consensus.cici.vcf -o tumor.gt.vcf

The output VCF has `AO` (SV-supporting reads) and `RO` (reference reads) per
breakpoint, read by SVCFit as `BPC` and `BEC` after scaling by mean read depth.

### 5. Germline heterozygous SNP detection (GATK4 + bcftools)

    gatk HaplotypeCaller -I normal.recal.bam -R ref.fa -O normal.germline.vcf.gz
    bcftools view -i 'GT="0/1"' normal.germline.vcf.gz -Oz -o normal.het.vcf.gz
    bcftools index -t normal.het.vcf.gz

### 6. Allele-specific copy-number profile (FACETS)

    snp-pileup -g -q15 -Q20 -P100 -r25,0 normal.het.vcf.gz tumor.snp.csv \
               normal.recal.bam tumor.recal.bam

``` r
library(facets)
rcmat <- readSnpMatrix("tumor.snp.csv")
xx    <- preProcSample(rcmat)
oo    <- procSample(xx, cval = 150)
fit   <- emcncf(oo)
write.table(fit$cncf, file = "tumor.facets.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

### 7. Per-SV breakpoint pileups at proximal heterozygous SNPs

    samtools view -b -f 1 -F 2 tumor.recal.bam --regions-file svs.bed > tumor.bp.bam
    samtools index tumor.bp.bam
    bcftools mpileup -f ref.fa -R svs_proximal_het.bed tumor.bp.bam > tumor.bp.pileup

### Caller-specific quirks

- **Manta**: writes `INFO/CIPOS` and `INFO/CIEND` natively — no reformatting
  needed.
- **Delly**: writes `INFO/CIPOS` but not `INFO/CIEND` — append `CIEND=-50,50`
  (or the per-call confidence interval if available).
- **GRIDSS**: uses paired-end / split-read counts in dedicated INFO fields;
  convert to LUMPY-style `INFO/MATEID`, `INFO/CIPOS`, and `INFO/CIEND` before
  SVtyper. A converter script `gridss_to_lumpy.py` is provided in `scripts/`.
- **Manta + Delly + GRIDSS consensus via SURVIVOR**: SURVIVOR drops
  `CIPOS`/`CIEND` from some merged records — run `add_cipos_ciend.py` to
  backfill them.

## Usage Workflow

The SVCFit pipeline consists of three main
steps: **Extraction**, **Characterization**, and **Calculation**.

### 1. Extract SV and SNP Information — `extract_info()`

Load and preprocess input VCF and CNV files. This step parses metadata,
processes breakends (BND), and handles heterozygous SNPs.
`extract_info()` internally performs:

1.  **Load input data** — `load_data()`
2.  **Process BND events** — `proc_bnd()`
3.  **Parse SV metadata** — `parse_sv_info()`
4.  **Parse heterozygous SNPs** — `parse_het_snp()`, `parse_snp_on_sv()`

``` r
info <- extract_info(
  p_het = "path/to/het_snps.vcf",
  p_onsv = "path/to/snps_on_sv.vcf",
  p_sv = "path/to/structural_variants.vcf",
  p_cnv = "path/to/cnv_file.txt",
  chr_lst = NULL,
  flank_del = 50, 
  QUAL_thresh = 100, 
  min_alt = 2, 
  tumor_only = FALSE
)
```

#### Function Arguments

#### `extract_info()`

| Argument | Type | Default | Description |
|----|----|----|----|
| `p_het` | Character | — | Path to VCF of heterozygous SNPs. |
| `p_onsv` | Character | — | Path to VCF of SNPs overlapping SV-supporting reads. |
| `p_sv` | Character | — | Path to SV VCF. |
| `p_cnv` | Character | — | Path to CNV file. |
| `chr_lst` | Character | NULL | Chromosomes to include. |
| `flank_del` | numeric | 50 | Max distance to consider deletion overlapping a BND. |
| `QUAL_thresh` | numeric | 100 | Minimum QUAL score. |
| `min_alt` | numeric | 2 | Minimum alternative reads. |
| `tumor_only` | Logical | FALSE | Whether SVs come from tumor-only calling. |

**Output:** A list of data frames containing parsed SV + SNP
information.

### 2. Annotate SVs Using CNV and SNP Information — `characterize_sv()`

This step integrates CNV and heterozygous SNPs to infer phasing,
zygosity, and overlapping CNV. `characterize_sv()` internally performs:

1.  **Assign SV IDs to SNPs** — `assign_svids()`
2.  **Summarizes phasing + zygosity** — `sum_sv_info()`
3.  **Assign CNV to SV** — `assign_cnv()`
4.  **Annotate overlapping CNV** — `annotate_cnv()`, `parse_snp_on_sv()`

``` r
sv_char <- characterize_sv(
  sv_phase = info$sv_phase, 
  sv_info = info$sv_info, 
  cnv = info$cnv,
  flank_snp = 500,
  flank_cnv = 1000
)
```

#### Function Arguments

#### `characterize_sv()`

| Argument    | Type       | Default | Description                       |
|-------------|------------|---------|-----------------------------------|
| `sv_phase`  | data.frame | —       | Phasing/zygosity from SNPs.       |
| `sv_info`   | data.frame | —       | Parsed SV metadata.               |
| `cnv`       | data.frame | —       | CNV data.                         |
| `flank_snp` | numeric    | 500     | Max assignment distance for SNPs. |
| `flank_cnv` | numeric    | 1000    | Max assignment distance for CNVs. |

### 3. Calculate SVCF for Structural Variants — `calculate_svcf()`

This step computes the **Structural Variant Cellular Fraction (SVCF)**.
and returns an annotated VCF file in data.frame format.

``` r
svcf_out <- calc_svcf(
  anno_sv_cnv = anno_sv_cnv,
  sv_info     = sv_info,
  thresh      = 0.1,
  samp        = "SampleID",
  exper       = "ExperimentID"
)
```

#### Function Arguments

| Argument      | Type       | Default | Description                            |
|---------------|------------|---------|----------------------------------------|
| `anno_sv_cnv` | data.frame | —       | CNV-annotated SVs.                     |
| `sv_info`     | data.frame | —       | Parsed SV info.                        |
| `thresh`      | numeric    | 0.1     | Threshold for SV-before-CNV inference. |
| `samp`        | character  | —       | Sample name.                           |
| `exper`       | character  | —       | Experiment name.                       |

**Output:** An annotated VCF-like data frame with additional fields for
VAF, Rbar, r, and SVCF.

1.  VAF: variant allele frequency
2.  Rbar: average break interval count in a sample
3.  r: inferred integer copy number of break intervals
4.  SVCF: structural variant cellular fraction.

### 4. Build tumor evolution tree — `build_tree()`

This step build the tumor evolutionary tree based on SV clusters
obtained from Dirichlet process Gaussian Mixture Model
(DP-GMM).Currently, this step is optimized for two sample longitudinal
data.

``` r
output <- cluster_data(
  pair_path  = pair_path,
  pur_path   = pur_path,
  data_dir   = data_dir,
  pair_num   = 1)
clones <- output[[3]]

build_tree(
  clones                    = clones,
  lineage_precedence_thresh = 0.2,
  sum_filter_thresh         = 0.2)
```

#### Function Arguments

cluster_data()

| Argument | Type | Default | Description |
|----|----|----|----|
| `pair_path` | character | — | Path to a tab-separated file with columns for ‘pre_BAT sample’ and ‘on_BAT sample’. |
| `pur_path` | character | \- | Path to a tab-separated file with columns for ‘sample’ and ‘purity’. |
| `data_dir` | character | — | Root directory containing SVCFit output BED files. |
| `pair_num` | numeric | 1 | The identifier (index or ID) for the specific sample pair (patient) being analyzed. |

build_tree()

| Argument | Type | Default | Description |
|----|----|----|----|
| `clones` | data.frame | — | SV clustering result. |
| `lineage_precedence_thresh` | numeric | 0.2 | Maximum violation of lineage precedence rule. |
| `sum_filter_thresh` | numeric | 0.2 | Maximum violation of sum condition rule |

**Output:** A tumor evolutionary tree rooted at the germline (G). Node
numbers correspond to SV cluster numbers. The branching depicts the
chronological occurrence of SV clusters.

### 5. Simulation & Benchmarking

SVCFit includes utility functions for processing simulation data from
VISOR and attaching “ground truth” labels to structural variants for
benchmarking.

5.1 read clonal assignment

``` r
truth <- load_truth(
  truth_path = "path/to/truth_beds", 
  overlap = FALSE
  )
```

This function has 1 arguments:

| Argument | Type | Default | Description |
|----|----|----|----|
| `truth_path` | Character | N/A | Path to BED files storing true structural variant information with clonal assignment. Each BED file should be named like `"c1.bed, c2.bed"`, etc for non-overlapping simulations and `"c11.bed, c22.bed"`, etc for overlapping simulations. Structural variants should be saved in separate BED files if they belong to different (sub)clones. |
| `overlap` | Logical | FALSE | Whether the simulation has SV-CNV overlap. |

The file path should follow this structure:

``` r
root/
├── true_clone/
│   ├── c1.bed/
│   ├── c2.bed/
│   ├── c3.bed/
│   └── .../
```

Parent nodes should always have lower number in name than its children
(i.e. c1.bed instead of c3.bed) and all child node bed file should
conatin its ancestors mutations.

<img src="inst/extdata/tree.png" width="487" />

5.2 attach clonal assignment to output

``` r
svcf_truth <- attach_truth(svcf_out, truth)
```

This function has 3 arguments:

| Variable | Type | Default | Description |
|----|----|----|----|
| `svcf_out` | DataFrame | N/A | The output from `calc_svcf` |
| `truth` | DataFrame | N/A | Stores the clone assignment for each structural variant designed in a simulation. |

This appends the known clonal assignment to the calculated SVCF output
for performance evaluation.

## Tutorial

``` r
library(SVCFit)
vignette("SVCFit_guide", package = "SVCFit")
```

## Frequently Asked Questions

**Can I use a different SV caller?**
Yes. Any caller that produces a per-SV breakpoint VCF compatible with SVtyper
will work. Differences across callers in breakpoint detection sensitivity,
split-read versus discordant-pair definitions, and quality filtering may produce
different REF/ALT counts and therefore different SVCF estimates from the same
data — see Discussion in the manuscript.

**Can I skip FACETS and use Battenberg / ASCAT instead?**
Yes, as long as you provide per-segment total copy number and gain/loss
classification in the same TSV format. A converter is in
`scripts/cnv_converter.R`.

**Do I need a matched normal?**
For the COMBAT analysis we used matched normals throughout. SVCFit can run on
tumor-only when a matched normal is unavailable, but tumor purity must be
supplied externally and ASCN inference becomes less reliable without
germline-heterozygous SNP calls; this is treated as an unsupported configuration
in the current release.

**What about complex SVs (chromothripsis, BFB, chromoplexy)?**
The closed-form SVCF estimators cover deletions, tandem duplications, inversions,
and the three classes of translocations. Multi-breakpoint complex SVs are not yet
explicitly modeled — the per-breakpoint estimates are still produced, but their
interpretation as a single cellular fraction is approximate. Future releases will
add structure-aware handling.

## Tool Versions

Tool versions used in the manuscript:

| Step | Tool | Version | Dataset |
|------|------|---------|---------|
| Read trimming | trim_galore | v0.6.1 | COMBAT |
| Alignment | bwa mem | v0.7.19 | All datasets |
| MarkDuplicates / BQSR | GATK4 | v4.6.2.0 | COMBAT |
| Somatic SV (single-caller) | Manta | v1.6.0 | All datasets |
| Somatic SV (multi-caller) | Manta + Delly v1.5.0 + GRIDSS v2.13.2 | — | COMBAT |
| Multi-caller merge | SURVIVOR | v1.0.7 | COMBAT |
| SV genotyping | SVtyper | v0.7.1 | All datasets |
| Germline SNP calling | GATK4 HaplotypeCaller | v4.6.2.0 | All datasets |
| Het-SNP filter | bcftools | v1.20 | All datasets |
| ASCN | FACETS | v0.6.2 | All datasets |
| Breakpoint read filter | samtools | v1.21 | All datasets |
| SNP pileup | bcftools mpileup | v1.20 | All datasets |

## Reference

1.  Cmero, Marek, Yuan, Ke, Ong, Cheng Soon, Schröder, Jan, Corcoran,
    Niall M., Papenfuss, Tony, et al., “Inferring Structural Variant
    Cancer Cell Fraction,” Nature Communications, 11(1) (2020), 730.
2.  Chen, X. et al. (2016) Manta: rapid detection of structural variants
    and indels for germline and cancer sequencing applications.
    Bioinformatics, 32, 1220-1222. <doi:10.1093/bioinformatics/btv710>
3.  Shen R, Seshan VE. FACETS: allele-specific copy number and clonal
    heterogeneity analysis tool for high-throughput DNA sequencing.
    Nucleic Acids Res. 2016 Sep 19;44(16):e131. doi: 10.1093/nar/gkw520.
    Epub 2016 Jun 7. PMID: 27270079; PMCID: PMC5027494.
4.  Auwera, Geraldine van der, and Brian D O’Connor. Genomics in the
    Cloud : Using Docker, GATK, and WDL in Terra. First edition.
    Sebastopol, CA: O’Reilly Media, 2020. Print.
5.  Li H. A statistical framework for SNP calling, mutation discovery,
    association mapping and population genetical parameter estimation
    from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93.
    doi: 10.1093/bioinformatics/btr509. Epub 2011 Sep 8. PMID: 21903627;
    PMCID: PMC3198575.
6.  Li, H., et al. (2009). The Sequence Alignment/Map format and
    SAMtools. *Bioinformatics*, Volume 25, Issue 16, August 2009, Pages
    2078–2079, <https://doi.org/10.1093/bioinformatics/btp352>
