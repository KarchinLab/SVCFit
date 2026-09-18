# SVCFit analysis workflow preservation

This directory holds the shared configuration framework for migrating the mixed legacy scripts currently in `SVCFit-2024-2026/script`. The original directory remains untouched. The file-by-file assessment, unavailable-source inventory, and migration recommendations are maintained with the resubmission working documents outside this repository.

To configure a machine:

```bash
cp config.example.sh config.local.sh
# Edit config.local.sh for the site and set SVCFIT_EXPECTED_COMMIT.
export SVCFIT_CONFIG=/absolute/path/to/config.local.sh
```

Preserved workflow scripts should source `lib/load_config.sh` and consume its exported variables. `config.local.sh`, logs, and outputs are ignored. No protected data, credentials, local paths, or generated results should be committed.

When approved workflow copies are migrated, replace their embedded paths with these variables and validate their outputs against the completed legacy runs on the protected server.
