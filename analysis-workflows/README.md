# SVCFit analysis workflow preservation

This directory holds the shared configuration and the maintained analysis-specific workflows migrated from the legacy script collection. The original directory remains outside this repository and is protected by the checksum manifest. The file-by-file assessment and correspondence are maintained with the resubmission working documents outside this repository.

To configure a machine:

```bash
cp config.example.sh config.local.sh
# Edit config.local.sh for the site and set SVCFIT_EXPECTED_COMMIT.
export SVCFIT_CONFIG=/absolute/path/to/config.local.sh
```

Each maintained workflow sources `lib/load_config.sh` and consumes its exported variables. `config.local.sh`, logs, and outputs are ignored. No protected data, credentials, local paths, or generated results should be committed.

## Migrated workflows

- `workflows/combat-germline-rerun`: PBMC-derived germline heterozygous-site filtering, tumor pileup at those sites, and paired SVCFit clustering/tree analysis.
- `workflows/combat-81875-downsampling`: 20-replicate 17x BAM downsampling, Delly/GRIDSS/Manta calling, SVTyper genotyping, SURVIVOR merging, and SVCFit-VCF preparation.

Both workflows pass local syntax, parse, configuration, and dry-run checks. They still require a one-sample smoke test and comparison with completed outputs on the protected server before being treated as production-validated.

The unresolved legacy `run_svcfit_boot.sh` was not installed as a maintained step because its comments and command disagree about whether the analysis is single-stage or longitudinal.
