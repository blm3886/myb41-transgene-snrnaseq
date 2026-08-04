# Cell Ranger workflow: methods log

Stage: align + quantify the 10x data against the custom reference

Last updated: 2026-07-20

---

## Input data (raw reads)

- **Location:** `/data2/nhartwick/law_lab_myb41_scrna/rawreads` (local on **seabiscuit**, 85 GB)
- **Platform:** 10x Genomics **3′ v3** — R1 = 28 bp (16 bp cell barcode + 12 bp UMI),
  R2 = 90 bp (cDNA), dual index (I1/I2 = 10 bp). All lane L003.
- **Samples (sample sheet):**| FASTQ prefix                                 | Sample name  | Role              |
  | -------------------------------------------- | ------------ | ----------------- |
  | `Col-0_S1_L003`                            | Col-0        | wild-type control |
  | `pFACT-MYB41_S2_L003`                      | pFACT-MYB41  | FACT construct    |
  | `pHORST-MYB41_S3_L003`                     | pHORST-MYB41 | HORST construct   |
  | Each sample has R1, R2, I1, I2`.fastq.gz`. |              |                   |

---

## Install Cell Ranger (one-time)

Cell Ranger is proprietary: downloaded from 10x Genomics.

- **Version:** cellranger-10.1.0
- **Downloaded:** signed URL from https://www.10xgenomics.com/support/software/cell-ranger/downloads
  (`curl -o cellranger-10.1.0.tar.gz "<signed-url>"`, 931 MB)
- **Installed to:** `/data1/bfernando/software/cellranger/cellranger-10.1.0/`
  ```bash
  mv cellranger-10.1.0.tar.gz /data1/bfernando/software/cellranger/
  cd /data1/bfernando/software/cellranger
  tar -xzvf cellranger-10.1.0.tar.gz
  ```
- **On PATH (persisted in ~/.bashrc):**
  ```bash
  export PATH=/data1/bfernando/software/cellranger/cellranger-10.1.0:$PATH
  ```
- Verify: `cellranger --version`  → `cellranger-10.1.0`

---

## Build the Cell Ranger reference (`mkref`)

Converts our combined FASTA + GTF (from notebook 02) into the indexed folder that

`cellranger count` requires. Builds a STAR index + gene lookups internally.

Inputs (from `reference/build/`):

- `Athaliana.Col-0.HPIv02.plus_pFusionMYB41.fasta` (8 seqs: Chr1-5, ChrC, ChrM, pFusionMYB41)
- `Athaliana.Col-0.HPIv02.plus_pFusionMYB41.gtf` (native + `AT4G28110.Fusion`)

```bash
cd /data1/bfernando/DE_project/reference/build
cellranger mkref \
  --genome=HPIv02_pFusionMYB41 \
  --fasta=Athaliana.Col-0.HPIv02.plus_pFusionMYB41.fasta \
  --genes=Athaliana.Col-0.HPIv02.plus_pFusionMYB41.gtf \
  --nthreads=16
```

- Output reference folder: `reference/build/HPIv02_pFusionMYB41/`
- Known risk: base annotation came via `gffread` (no `gene` lines, uses `transcript`)
- while the contig GTF (AGAT) has `gene`/`mRNA`/nonstandard `insertion`.
- If mkref rejects it, harmonize the GTF and re-run.
  - STATUS: **PASSED (2026-07-21)** after a fix. First run FAILED with
    `GexReferenceError: Transcripts with no exons were detected`.
  - Root cause: the source Araport11 GFF3 has **no `exon` features** (only gene/mRNA/CDS/UTR), and `gffread -T` did not synthesize exons for ~3,824 CDS-only (UTR-less) transcripts.
  - FIX: re-convert with **`gffread --force-exons`** (0 transcripts without exons), re-merge with contig
    GTF, `rm -rf mkref_HPIv02_pFusionMYB41`, re-run.
  - Notebook 02 Step 3 patched to include `--force-exons`. Reference now at `reference/build/HPIv02_pFusionMYB41/`.

---

## Quantify each sample (`cellranger count`)

Run once per sample (Col-0, pFACT, pHORST) against the `mkref` reference.

```bash
cellranger count \
  --id=<sample_id> \
  --transcriptome=/data1/bfernando/DE_project/reference/build/HPIv02_pFusionMYB41 \
  --fastqs=/data2/nhartwick/law_lab_myb41_scrna/rawreads \
  --sample=<fastq_prefix e.g. Col-0 | pFACT-MYB41 | pHORST-MYB41> \
  --chemistry=auto \
  --include-introns=true \
  --create-bam=true \
  --localcores=16 --localmem=64
```

- **⚠ CHEMISTRY MISMATCH (2026-07-21):** first Col-0 run used `--chemistry=SC3Pv3` and
  FAILED — **Estimated cells = 54, Valid Barcodes = 3.4%**.
- Root cause: the data is actually **3' v4 / GEM-X** (`SC3Pv4`), not v3.
- Confirmed by matching R1 barcodes vs Cell Ranger whitelists: 78% match to `3M-3pgex-may-2023` (v4) vs 0.5% to v3.
  v3 and v4 both have R1 = 28 bp, so read length can't distinguish them — only the
  barcode whitelist does.
- **FIX: use `--chemistry=auto`** (detects SC3Pv4).
- Delete the failed `counts/Col-0` before re-running
