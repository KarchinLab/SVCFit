# COMBAT sample-81875 downsampling workflow

This workflow downsamples sample 81875 to 17x across 20 seeded replicates, calls SVs with Delly, GRIDSS, and Manta, genotypes the calls with SVTyper, merges calls with SURVIVOR, and prepares one SVCFit-compatible VCF per replicate.

The maintained downsampling step is based on the safer legacy implementation formerly named `old_down_cov.sh`. It validates the sampling fraction, records the seed, replaces read groups, checks the final BAM, and optionally removes intermediates.

## Run

Copy `../../config.example.sh` to `../../config.local.sh`, edit the protected-server values, and export `SVCFIT_CONFIG`.

```bash
./submit.sh --single --dry-run  # inspect a one-replicate dependency chain
./submit.sh --single            # protected-server smoke test
./submit.sh                     # configured number of replicates
```

The submission wrapper creates the workflow log directory and uses `afterok` dependencies. Set `DOWNSAMPLE_REPLICATES`, sample names, coverage values, seeds, tool locations, and recovery inputs in the central configuration.

The configured environments must provide `samtools`, `delly`, `bcftools`, GRIDSS, Manta, `svtyper-sso`, SURVIVOR, Python, and R with `optparse`, `dplyr`, `tidyr`, and `stringr`. Caller-specific VCF layouts are inherited from the recorded legacy toolchain and must be confirmed during the protected-server smoke test.

The pipeline ends with SVCFit-ready VCF preparation. The legacy `run_svcfit_boot.sh` is intentionally excluded because its documentation says `--stages svcfit` and same-sample inputs while its command runs the default stages with a longitudinal partner. That intent must be resolved before adding a downstream inference step.

The scripts have been checked locally for syntax, R parsing, configuration loading, and submission dry runs. Caller execution and scientific comparison with legacy outputs remain protected-server validation tasks.
