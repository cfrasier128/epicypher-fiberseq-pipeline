# Fiberseq Nextflow Pipeline

Summary

- Purpose: Align PacBio BAMs, optionally run variant calling + phasing,
  add MSP/nucleosome annotations, run fiberseq QC, and optionally create
  pileups + BigWig tracks.
- Entry point: `main.nf` (DSL2 workflow) plus `subworkflows/fiberseq-qc.nf`.
- Convenience: `Makefile` provides install/lint/format/clean and
  reference-prep helpers.

Requirements

- Nextflow (the `Makefile` can download a pinned version).
- Containers:
  - `-profile local` uses Docker.
  - `-profile slurm` uses Singularity/Apptainer.
- `samtools` is required by the reference prep script.

Inputs

Sample sheet (required)

- `--sample_sheet`: TSV with header and columns:
  `samp_name <TAB> bam_path <TAB> ref_name`.
- Multiple rows can share the same `samp_name` (those BAMs will be aligned
  individually then merged).

Example: `inputs/sample_sheet.tsv`

```tsv
samp_name\tbam_path\tref_name
sampleA\tinputs/bams/a1/sampleA.bam\thg38
sampleA\tinputs/bams/a2/sampleA.bam\thg38
sampleB\tinputs/bams/b1/sampleB.bam\tchm13
```

Reference sheet (required)

- `--ref_sheet_path`: TSV with header and columns:
  `ref_name <TAB> ref_fasta <TAB> ref_index`.
- The `ref_name` values must match between the sample sheet and reference sheet.

Example `inputs/reference_sheet.tsv`:

```tsv
ref_name\tref_fasta\tref_index
hg38\t/refs/hg38/hg38.fasta\t/refs/hg38/hg38.fasta.fai
chm13\t/refs/chm13/chm13.fasta\t/refs/chm13/chm13.fasta.fai
```

Parameters (defaults shown)

- `--outdir` (default `${workflow.launchDir}/results`) — top-level output
  directory.
- `--confidence_ml_val` (default `250`) — ML threshold for `ft add-nucleosomes`
  and pileups.
- `--minimum_msp_dist` (default `10`) — MSP length filter used for pileups
  (`ftx "len(msp) > ..."`).

Optional steps

- `--pb_qc` (default `false`) — generate PacBio QC reports.
- `--phase_reads` (default `false`) — run `deepvariant` + `sawfish` then
  `hiphase` haplotagging.
- `--create_bigwigs` (default `false`) — create pileup TSVs and BigWigs.
- `--debug` (default `false`) — prints helpful channel `view()` messages.

Makefile shortcuts

- Install Nextflow locally into the repo:

```bash
make install
```

- Create a blank sample sheet template:

```bash
make sample_sheet_template
```

- Prepare a reference FASTA + index and append to a reference sheet:

```bash
make references ref=hg38
make references ref=chm13
```

This writes FASTA files under `references/<ref>/` and creates/appends to
`inputs/reference_sheet.tsv` (configurable via `REFDIR` and `REFSHEET`).

- Lint / auto-format the Nextflow scripts:

```bash
make check
make format
```

High-level workflow steps (what `main.nf` does)

- Read the sample sheet and the reference sheet and join by `ref_name`.
- `align_bams`: align each input BAM with `pbmm2 align` (HiFi preset),
  producing `*.aligned.bam` + `*.bai`.
- `merge_bams`: group by `samp_name` and merge aligned BAMs into
  `${samp_name}.${ref_name}.aligned.bam` (then `samtools index`).
- If `--pb_qc`: `pacbio_qc` writes reports under
  `${outdir}/sequencing_qc/${samp_name}/`.
- If `--phase_reads`:
  - `deepvariant` produces `${samp_name}.${ref_name}.deepvariant.vcf.gz`.
  - `sawfish` produces `${samp_name}.${ref_name}.structural_variants.vcf.gz`.
  - `hiphase` produces `${samp_name}.${ref_name}.haplotagged.bam` and phased
    VCFs under `${outdir}/phased_output/${samp_name}/`.
- `call_msps`: runs `ft add-nucleosomes` and publishes annotated BAMs under
  `${outdir}/fire_bams/${samp_name}/`.
- `fiberseq_qc_workflow`: runs Stergachis-style fiberseq QC and publishes under
  `${outdir}/fiberseq-qc/${sample_id}/`.
- If `--create_bigwigs`: creates `${samp_name}.pileup_all.tsv.gz` and converts
  to BigWigs under `${outdir}/pileups/${samp_name}/`.

Published output locations

- `${outdir}/sequencing_qc/<sample>/` — PacBio QC PDFs/plots
  (only if `--pb_qc`).
- `${outdir}/phased_output/<sample>/` — haplotagged BAM + phased VCFs
  (only if `--phase_reads`).
- `${outdir}/fire_bams/<sample>/` — `*.6ma.nucs.bam` and index.
- `${outdir}/fiberseq-qc/<sample>/` — QC tables + PDFs.
- `${outdir}/pileups/<sample>/` — `*.pileup_all.tsv.gz` and `*.bw` tracks
  (only if `--create_bigwigs`).

Example runs

Local run (Docker), align + MSP/QC (minimum useful run)

```bash
nextflow run main.nf \
  --sample_sheet inputs/sample_sheet.tsv \
  --ref_sheet_path inputs/reference_sheet.tsv \
  -profile local
```

Local run with phasing + bigwigs

```bash
nextflow run main.nf \
  --sample_sheet inputs/sample_sheet.tsv \
  --ref_sheet_path inputs/reference_sheet.tsv \
  --pb_qc true \
  --phase_reads true \
  --create_bigwigs true \
  --outdir results/full_run \
  -profile local
```

SLURM run (Singularity)

```bash
nextflow run main.nf \
  --sample_sheet inputs/sample_sheet.tsv \
  --ref_sheet_path inputs/reference_sheet.tsv \
  -profile slurm
```

Tips

- Use absolute paths in `inputs/reference_sheet.tsv` (the
  `prepare_references.sh`
  script writes canonical paths for this reason).
- If you enable `-profile debug`, Nextflow will emit
  trace/timeline/report/dag files; use `--debug true` if you also want the
  channel `view()` messages.
