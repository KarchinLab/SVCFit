# Bootstrap SV Pipeline — Session Notes

## Overview

Pipeline for sample **81875**: downsample a tumor BAM to 17x coverage across 20 bootstrap replicates, call SVs with three callers, genotype with SVtyper, merge with SURVIVOR, and prepare output for SVCFit analysis.

---

## Script Directory

```
/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/script/sv_call/
```

Scripts present:
```
del_svtyp.sh            run_down_samp_gridss.sh
down_cov_81875.sh       run_down_samp_manta.sh
get_read_for_surv.R     run_pipeline.sh
gri_svtyp.sh            run_survivor.sh
man_svtyp.sh            run_svcfit_cli.R
proc_surv.sh            surv2svcfit.R
run_down_samp_delly.sh  svcfit.sh
```

---

## Pipeline Dependency Chain

```
downsample (array 1-20)
  ├── delly  (array 1-20) → del_svtyp (array 1-20) ──┐
  ├── gridss (array 1-20) → gri_svtyp (array 1-20) ──┼── survivor (array 1-20) → proc_surv (array 1-20)
  └── manta  (array 1-20) → man_svtyp (array 1-20) ──┘
```

SLURM chaining uses `--dependency=afterok`: if any replicate in an array step fails, all downstream jobs are held (correct behavior — avoids running SURVIVOR on incomplete data).

---

## Changes Made to Scripts

### All scripts: array range updated to `1-20`

| Script | Before | After |
|---|---|---|
| `run_down_samp_delly.sh` | `--array=1-30%10` | `--array=1-20%10` |
| `run_down_samp_gridss.sh` | `--array=6-30%10` | `--array=1-20%10` |
| `run_down_samp_manta.sh` | `--array=6-30%10` | `--array=1-20%10` |
| `del_svtyp.sh` | `--array=1-30` | `--array=1-20` |
| `gri_svtyp.sh` | `--array=6-30` | `--array=1-20` |
| `man_svtyp.sh` | `--array=6-30` | `--array=1-20` |
| `run_survivor.sh` | `--array=1-30` | `--array=1-20` |
| `proc_surv.sh` | `--array=1-30` | `--array=1-20` |

### All downstream scripts: BAM filename corrected

All 7 downstream scripts had `81875_puritymatched_17x_boot$id.rg.bam` replaced with `81875_covdown_17x_boot$id.rg.bam` to match the actual output of `down_cov_81875.sh`.

### Log paths: relative → absolute

Five scripts used relative `log/` paths in their `#SBATCH --output/--error` directives, which would resolve unpredictably depending on the submission directory. Changed to absolute paths:

| Script | Log path |
|---|---|
| `down_cov_81875.sh` | `/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/script/sv_call/log/` |
| `del_svtyp.sh` | same |
| `gri_svtyp.sh` | same |
| `man_svtyp.sh` | same |
| `proc_surv.sh` | same |

### `proc_surv.sh`: R script path corrected

```bash
# Before
Rscript $sv_dir/survivor/surv2svcfit.R ...
# ($sv_dir = .../svcfit/sv_call — wrong location)

# After
Rscript /projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/script/sv_call/surv2svcfit.R ...
```

---

## New File: `run_pipeline.sh`

Wrapper that submits the full pipeline as a chained set of SLURM array jobs.

### Usage

```bash
bash run_pipeline.sh                     # submit all 20 replicates
bash run_pipeline.sh --single            # submit replicate 1 only (for testing)
bash run_pipeline.sh --dry-run           # print sbatch commands without submitting
bash run_pipeline.sh --single --dry-run  # preview single-replicate chain
```

### How it works

- Uses `sbatch --parsable` to capture each job ID.
- Passes the ID to the next step via `--dependency=afterok:<JID>`.
- SURVIVOR waits for all three SVtyper jobs: `--dependency=afterok:<DEL>:<GRI>:<MAN>`.
- `--single` injects `--array=1` into every `sbatch` call, overriding the script's own `#SBATCH --array`.
- `--dry-run` replaces every `sbatch` call with an `echo` and assigns sequential fake job IDs so the dependency chain can be inspected without touching the cluster.
- Creates all required log directories before any job is submitted.

### Dry-run output (representative)

```
sbatch --parsable down_cov_81875.sh                                         → JID 1001
sbatch --parsable --dependency=afterok:1001 run_down_samp_delly.sh          → JID 1002
sbatch --parsable --dependency=afterok:1001 run_down_samp_gridss.sh         → JID 1003
sbatch --parsable --dependency=afterok:1001 run_down_samp_manta.sh          → JID 1004
sbatch --parsable --dependency=afterok:1002 del_svtyp.sh                    → JID 1005
sbatch --parsable --dependency=afterok:1003 gri_svtyp.sh                    → JID 1006
sbatch --parsable --dependency=afterok:1004 man_svtyp.sh                    → JID 1007
sbatch --parsable --dependency=afterok:1005:1006:1007 run_survivor.sh       → JID 1008
sbatch --parsable --dependency=afterok:1008 proc_surv.sh                    → JID 1009
```

---

## Required Inputs (pre-existing)

### Data files

| Variable | Path |
|---|---|
| Tumor BAM | `/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/bam_hg38/81875_hg38_recal_sorted.bam` (+ `.bai`) |
| Normal BAM | `/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/bam_hg38/PBMC_10_hg38_recal_sorted.bam` (+ `.bai`) |
| Reference genome | `/projects/karchin-lab-hpc/hg38_reference/hg38.analysisSet.fa` (+ `.fai`, `.dict`) |
| Reference Delly VCF *(optional, for SV recovery)* | `/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/sv_call/delly/d81875/81875.vcf` |
| Recovery file *(optional, per replicate)* | `/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/sv_call/survivor/recover_sv/boot{1..20}` |

The recovery file and reference Delly VCF are only used if the recovery file exists and is non-empty (`surv2svcfit.R` checks `has_nonempty_file()` and skips the block otherwise).

### Tool paths

| Tool | Path / expectation |
|---|---|
| `samtools` | in `$PATH` |
| `delly` | in `$PATH` |
| `bcftools` | in `$PATH` |
| `gridss` | in `$PATH` |
| `GRIDSS jar` | `/home/yliu498/miniforge3/envs/gridss/share/gridss-2.13.2-6/gridss.jar` |
| Manta | `~/manta/bin/configManta.py` and `~/manta/libexec/convertInversion.py` |
| `samtools` (conda) | `~/.conda/envs/visor/bin/samtools` (used by Manta's convertInversion) |
| `svtyper-sso` | in `$PATH` |
| `SURVIVOR` | in `$PATH` |
| `Rscript` | in `$PATH` |

### Log directories (created automatically by `run_pipeline.sh`)

```
/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/script/sv_call/log/
/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/sv_call/delly/log/
/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/sv_call/gridss/log/
/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/sv_call/manta/log/
/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/svcfit/sv_call/survivor/log/
```

---

## Data Flow (output → input consistency check)

All handoffs verified consistent. Summary:

| Step | Produces | Consumed by |
|---|---|---|
| `down_cov_81875.sh` | `.../down_samp/bootstrap_$id/81875_covdown_17x_boot$id.rg.bam` | delly, gridss, manta, all svtyper scripts |
| `run_down_samp_delly.sh` | `.../delly/boot$id/p_boot$id.vcf` | `del_svtyp.sh` |
| `run_down_samp_gridss.sh` | `.../gridss/boot$id/p_boot$id.vcf` | `gri_svtyp.sh` |
| `run_down_samp_manta.sh` | `.../manta/boot$id/results/variants/p_boot$id.vcf` | `man_svtyp.sh` |
| `del_svtyp.sh` | `.../svtyp/boot$id/delly_boot$id.vcf` | `run_survivor.sh`, `proc_surv.sh` |
| `gri_svtyp.sh` | `.../svtyp/boot$id/gridss_boot$id.vcf` | `run_survivor.sh`, `proc_surv.sh` |
| `man_svtyp.sh` | `.../svtyp/boot$id/manta_boot$id.vcf` | `run_survivor.sh`, `proc_surv.sh` |
| `run_survivor.sh` | `.../survivor/boot$id/boot$id_wbnd.vcf` | `proc_surv.sh` |
| `proc_surv.sh` | `.../survivor/boot$id/svcfit_boot$id.vcf` | SVCFit analysis |

---

## Testing Strategy

### Level 1 — Syntax check (instant, no cluster needed)

```bash
for f in down_cov_81875.sh run_down_samp_delly.sh run_down_samp_gridss.sh \
          run_down_samp_manta.sh del_svtyp.sh gri_svtyp.sh man_svtyp.sh \
          run_survivor.sh proc_surv.sh run_pipeline.sh; do
    echo -n "Checking $f ... "
    bash -n $f && echo "OK"
done
```

### Level 2 — Dry-run (validates submission logic, no compute)

```bash
bash run_pipeline.sh --dry-run           # all 20 replicates
bash run_pipeline.sh --single --dry-run  # replicate 1 only
```

SLURM's own validator (checks resource limits and queue):
```bash
sbatch --test-only down_cov_81875.sh
```

### Level 3 — Single replicate functional test

```bash
bash run_pipeline.sh --single
```

Monitor with:
```bash
squeue -u $USER
```

Once replicate 1 completes cleanly, run the full pipeline:
```bash
bash run_pipeline.sh
```
