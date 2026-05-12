# SVCFit 1.0.0

* Initial release.
* Estimates Structural Variant Cellular Fraction (SVCF) for inversions,
  deletions, tandem duplications, and translocations.
* Integrates SV calls with allele-specific copy number variation (FACETS)
  and heterozygous germline SNPs for phasing and zygosity inference.
* Dirichlet-process Gaussian mixture model (DP-GMM) clustering of SV
  cellular fractions across paired longitudinal samples.
* Phylogenetic spanning-tree inference via a modified Gabow-Myers algorithm
  with CCF sum-condition filtering and linear-chain penalty.
* Utility functions for VISOR simulation benchmarking (`load_truth`,
  `attach_truth`) and SVclone comparison (`analyze_svclone`).
