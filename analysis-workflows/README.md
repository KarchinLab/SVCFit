# SVCFit analysis workflow preservation

This directory holds the shared configuration and migration design for the mixed legacy scripts currently in `SVCFit-2024-2026/script`. The original directory remains untouched. See `SCRIPT_DIRECTORY_REVIEW.md` for the file-by-file assessment, unavailable-source inventory, recommended repository layout, and migration order.

To configure a machine:

```bash
cp config.example.sh config.local.sh
# Edit config.local.sh for the site and set SVCFIT_EXPECTED_COMMIT.
export SVCFIT_CONFIG=/absolute/path/to/config.local.sh
```

Preserved workflow scripts should source `lib/load_config.sh` and consume its exported variables. `config.local.sh`, logs, and outputs are ignored. No protected data, credentials, local paths, or generated results should be committed.

The configuration and review are the first preservation commit. Workflow copies should be migrated in the order described in `SCRIPT_DIRECTORY_REVIEW.md`, with their embedded paths replaced by these variables and their outputs validated against the completed legacy runs on the protected server.
