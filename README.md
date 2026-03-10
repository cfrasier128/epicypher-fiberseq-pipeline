# Fiber-seq Nextflow Pipeline
<h3>Disclaimer</h3>
While this repository is maintained by EpiCypher, Inc., we make no promises to troubleshoot or offer technical support. Please report any bugs encountered, however we make no promises to fix them.

---

<h3>Summary</h3>

- Purpose: This is a Nextflow based analysis pipeline focused on analyzing Fiber-seq data from a PacBio sequencing instrument. This pipeline will carry unaligned bams through alignment, merging, fibertools, and fiberseq-qc. Optionally, FIRE peak calling, haplotype phasing, and bigwig creation are supported. 
- Entry point: `main.nf` (DSL2 workflow).
- Convenience: `Makefile` provides install/lint/format/clean and
  reference-prep helpers.

---

<h3>Requirements</h3>

- Nextflow (the `Makefile` can download a pinned version). It is also possible to install via Conda.
- Containers:
  - `-profile local` uses Docker.
  - `-profile slurm` uses Singularity/Apptainer. (Can also be used to support Docker)
- `samtools` is required by the reference prep script.

---

<h3>Inputs</h3>

<h4>Required</h4>

- `--sample_sheet`: TSV with header and columns:
  `samp_name <TAB> bam_path <TAB> ref_name`.
  - Multiple rows can share the same `samp_name` (those BAMs will be aligned individually then merged, useful for technical sequencing replicates).\


- `--ref_sheet_path`: TSV with header and columns: `ref_name <TAB> fasta_path <TAB> fasta_index`.
  - Can be created using `prepare_references.sh`
  - The `ref_name` values must match between the sample sheet and reference sheet.

- `-profile`: Nextflow explicit parameter, determines method for job execution. Use one of the following:
  - local: executes job using local resources and Docker as container method
  - aws_env: executes job using Amazon Web Service Batch compute environment
  - slurm: executes job locally using slurm as job scheduler. Singularity or Docker can be used as container method.


Example: `inputs/sample_sheet.tsv`

| samp_name | bam_path | ref_name |
|-----------|----------|----------|
| sampleA | inputs/bams/a1/sampleA.bam | hg38 |
| sampleA | inputs/bams/a2/sampleA.bam | hg38 |
| sampleB | inputs/bams/b1/sampleB.bam | chm13 |


Example `inputs/reference_sheet.tsv`:

| ref_name | ref_fasta | ref_index |
|----------|-----------|-----------|
| hg38 | /refs/hg38/hg38.fasta | /refs/hg38/hg38.fasta.fai |
| chm13 | /refs/chm13/chm13.fasta | /refs/chm13/chm13.fasta.fai |


<h4>Optional Parameters</h4>

- `--outdir` (default `${workflow.launchDir}/results`) — top-level output
  directory.
- `--confidence_ml_val` (default `250`) — ML threshold to use for both 6mA and 5mC for `ft add-nucleosomes`
  and pileups.
- `--minimum_msp_dist` (default `10`) — MSP length filter used for pileups
  (`ftx "len(msp) > ..."`).

<h4>Optional steps</h4>

- `--pb_qc` (default `false`) — generate PacBio QC reports.
- `--phase_reads` (default `false`) — run `deepvariant` + `sawfish` then
  `hiphase` haplotagging.
- `--create_bigwigs` (default `false`) — create pileup TSVs and BigWigs.
- `--debug` (default `false`) — prints helpful channel `view()` messages.

---

<h3>Makefile shortcuts</h3>

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
`inputs/reference_sheet.tsv` (configurable via `REFDIR` and `REFSHEET` in the Makefile).

- Lint / auto-format the Nextflow scripts:

```bash
make check
make format
```

---

<h3>High-level workflow steps</h3>

Short descriptions of each Nextflow job step:
- Read the sample sheet and the reference sheet and join by `ref_name`.
- `align_bams`: align each input BAM with `pbmm2 align` (HiFi preset).
- `merge_bams`: group by `samp_name` from first column of sample sheet and merge aligned BAMs into (then index with `samtools index`).
- If using `--pb_qc`: creates sample level sequencing QC reports.
- If using `--phase_reads`:
  - `deepvariant` produces short nucleotide polymorphism (SNP/SNV) calls.
  - `sawfish` produces structural variant (SV) calls.
  - `hiphase` uses both variant calls (SNP and SV) to haplotype phase reads.
- `call_msps`: Paint bams with MSP and nucleosome calls.
- `fiberseq_qc_workflow`: runs Stergachis-style fiberseq QC.
- If using `--create_bigwigs`: creates methylation (6mA/5mC) and nucleosome pileups and converts them to BigWigs.

---

<h3>Published output locations</h3>

- `${outdir}/1_fire_bams/<sample>/` — Fibertools labelled bam and index.
- `${outdir}/2_fiberseq-qc/<sample>/` — Fiber-seq QC tables + PDFs.
- `${outdir}/3_phased_output/<sample>/` — haplotagged BAM + phased VCFs (only if `--phase_reads`).
- `${outdir}/4_pileups/<sample>/` — Pileups and bigwigs tracks (only if `--create_bigwigs`).
- `${outdir}/5_sequencing_qc/<sample>/` — Per sample PacBio QC PDFs/plots (only if `--pb_qc`).
- `${outdir}/6_FIRE_peaks/<sample>/` — Per sample PacBio QC PDFs/plots (only if `--peak_call`).

---


<h3>Example runs</h3>

Local run (Docker), align + MSP/QC (minimum required parameters)

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

---

Tips

- Use absolute paths in `inputs/reference_sheet.tsv` (the
  `prepare_references.sh`
  script writes canonical paths for this reason).
- If you enable `-profile debug`, Nextflow will emit
  trace/timeline/report/dag files; use `--debug true` if you also want the
  channel `view()` messages.

