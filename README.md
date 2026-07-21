# MYB41 transgene snRNA-seq

Single-nucleus RNA-seq (10x Genomics) to detect **promoter::MYB41** transgene expression
in *Arabidopsis* (Col-0), using a custom reference (Col-0 genome + a `pFusionMYB41`
transgene contig) and detecting a transcribed cloning scar unique to the transgene.

- **Constructs:** `pFACT::MYB41` and `pHORST::MYB41` (+ Col-0 wild-type control)
- **Platform:** 10x Chromium **3′ v4 (GEM-X)** — R1 = 28 bp (16 bp barcode + 12 bp UMI), R2 = 90 bp
- **Detection marker:** "Scar 2", a transcribed cloning scar at contig positions 390–417

## Repository contents

| Path | What it is |
|---|---|
| `notebooks/01_compare_constructs_2.ipynb` | Prove FACT vs HORST differ only in the promoter; validate Scar 2 |
| `notebooks/02_build_reference.ipynb` | Build the custom Col-0 HPIv02 + pFusionMYB41 reference (FASTA + GTF) |
| `notebooks/03_raw_reads_qc.ipynb` | Raw-FASTQ QC (quality, chemistry, contamination) |
| `results/results1_construct_observations.md` | Written summary of the construct findings |
| `results/qc/…web_summary.html` | Cell Ranger QC report(s) |
| `notes/03_cellranger_workflow.md` | Detailed methods log |

> Large data (genome, reference index, FASTQs, BAMs, count matrices) is **not** in the repo
> (see `.gitignore`); it lives on the analysis server.

---

## Environment & tools

- **Cell Ranger** `10.1.0` (installed at `/data1/bfernando/software/cellranger/cellranger-10.1.0`)
- **conda env `de_project`** — python 3.10, biopython, pandas, numpy, matplotlib, pysam,
  samtools 1.24, gffread 0.12.9 (for reference prep)

---

## Pipeline & exact commands

### Step 1 — Build the combined reference FASTA + GTF (notebook 02)

Assembled in `notebooks/02_build_reference.ipynb`:
1. Fix the malformed contig FASTA (add a `>pFusionMYB41` header).
2. Append the contig to the Col-0 **HPIv02** genome → combined FASTA (8 sequences:
   Chr1–5, ChrC, ChrM, pFusionMYB41).
3. Convert the base annotation **GFF3 → GTF** and merge with the contig GTF:

```bash
# the source GFF3 has NO exon features (only gene/mRNA/CDS/UTR); --force-exons is REQUIRED,
# otherwise ~3,824 CDS-only transcripts get no exon line and Cell Ranger mkref rejects the GTF.
gffread Athaliana.Col-0.HPIv02.gene.gff3 -T --force-exons -o Athaliana.Col-0.HPIv02.gene.gtf
cat Athaliana.Col-0.HPIv02.gene.gtf Athaliana.Col-0.HPIv01.pFusionMYB41_only.gtf \
    > Athaliana.Col-0.HPIv02.plus_pFusionMYB41.gtf
```

### Step 2 — Build the Cell Ranger reference index (`mkref`)

```bash
cd /data1/bfernando/DE_project/reference/build
cellranger mkref \
  --genome=HPIv02_pFusionMYB41 \
  --fasta=Athaliana.Col-0.HPIv02.plus_pFusionMYB41.fasta \
  --genes=Athaliana.Col-0.HPIv02.plus_pFusionMYB41.gtf \
  --nthreads=16
```

| Flag | Value | Purpose |
|---|---|---|
| `--genome` | `HPIv02_pFusionMYB41` | name of the output reference folder |
| `--fasta` | combined genome FASTA | Col-0 chromosomes + transgene contig (8 sequences) |
| `--genes` | combined GTF | native + `AT4G28110.Fusion` gene models |
| `--nthreads` | `16` | CPUs for STAR index generation |

Output: `reference/build/HPIv02_pFusionMYB41/` (contains `reference.json`, `star/`, `genes/`, `fasta/`).

### Step 3 — Generate the count matrix (`cellranger count`)

Run **once per sample** (Col-0, pFACT-MYB41, pHORST-MYB41), changing `--id` and `--sample`:

```bash
cd /data1/bfernando/DE_project/counts
cellranger count \
  --id=Col-0_v4 \
  --transcriptome=/data1/bfernando/DE_project/reference/build/HPIv02_pFusionMYB41 \
  --fastqs=/data2/nhartwick/law_lab_myb41_scrna/rawreads \
  --sample=Col-0 \
  --chemistry=auto \
  --include-introns=true \
  --create-bam=true \
  --localcores=16 --localmem=64
```

| Flag | Value | Purpose |
|---|---|---|
| `--id` | `Col-0_v4` | unique run id / output folder |
| `--transcriptome` | `…/HPIv02_pFusionMYB41` | the `mkref` reference from Step 2 |
| `--fastqs` | `/data2/nhartwick/law_lab_myb41_scrna/rawreads` | folder of FASTQs |
| `--sample` | `Col-0` \| `pFACT-MYB41` \| `pHORST-MYB41` | FASTQ filename prefix (before `_S#`) |
| `--chemistry` | `auto` | **auto-detect** — correctly identifies 3′ v4 (see note below) |
| `--include-introns` | `true` | count intronic reads (default true; required for single-nucleus) |
| `--create-bam` | `true` | emit the aligned BAM (needed to inspect Scar-2 reads) |
| `--localcores` / `--localmem` | `16` / `64` | resource caps (shared server) |

Output per sample: `<id>/outs/` with `web_summary.html`, `filtered_feature_bc_matrix/`,
`metrics_summary.csv`, and `possorted_genome_bam.bam`.

---

## Important notes / gotchas

- **Chemistry is 3′ v4, not v3.** v3 and v4 share the same read structure (R1 = 28 bp), so
  read length can't distinguish them — only the barcode whitelist can. Forcing `--chemistry=SC3Pv3`
  **fails** (valid barcodes 3.4%, 54 cells). Use `--chemistry=auto` (detects `SC3Pv4`).
- **`--force-exons` is required** when converting the source GFF3 (it has no `exon` features);
  without it, Cell Ranger `mkref` errors with "transcripts with no exons".
- **Reference validated:** in the Col-0 WT control, `AT4G28110.Fusion` = 0 counts (correct
  negative control) and 20,547 genes are detected with a sensible expression profile.
