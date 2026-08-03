**Cell Ranger :**

```bash
/data1/bfernando/software/cellranger
```

1. Custom Reference Building : /data1/bfernando/DE_project/notebooks/02_build_reference.ipynb
   builds the combined FASTA (Col-0 HPIv02 genome + `pFusionMYB41` transgene contig) and the
   merged GTF. Then `mkref` turns those into the indexed folder `cellranger count` needs:

   ```bash
   cellranger mkref --genome=HPIv02_pFusionMYB41 \
     --fasta=Athaliana.Col-0.HPIv02.plus_pFusionMYB41.fasta \
     --genes=Athaliana.Col-0.HPIv02.plus_pFusionMYB41.gtf --nthreads=16
   ```

   Output: `reference/build/HPIv02_pFusionMYB41/` — this is what `--transcriptome` points at below.
2. Aligned Reads to Ref Genome and generated counts file :

   1. Counts File at : `SAMPLE_OUTS = "/data1/bfernando/DE_project/Col-0_v4/outs`
   2. ```
      cellranger count --id=Col-0_v4 \
        --transcriptome=/data1/bfernando/DE_project/reference/build/HPIv02_pFusionMYB41 \
        --fastqs=/data2/nhartwick/law_lab_myb41_scrna/rawreads \
        --sample=Col-0 \
        --chemistry=auto \
        --include-introns=true \
        --create-bam=true --localcores=16 --localmem=64
      ```

   `--sample=Col-0` --> matches the sample-name field in the Cell Ranger filename convention
   `{sample}_S{n}_L{lane}_{read}_001.fastq.gz`, so only the 4 `Col-0_*` files were used
   (WT only). The 8 `pFACT-MYB41_*` / `pHORST-MYB41_*` files were ignored.

   `--chemistry=auto` --> because the first run used `--chemistry=SC3Pv3` and failed with only
   **54 cells** (1.2% valid barcodes). Diagnosed by barcode-whitelist matching:
   **78%** matched `3M-3pgex-may-2023` (v4) vs **0.5%** to the v3 list.
   Both v3 and v4 have R1 = 28 bp, so read length cannot distinguish them — only the
   barcode whitelist can. `auto` lets Cell Ranger detect SC3Pv4 from the barcodes.
3. The run detected **864 cells**. The barcode-rank plot
   (/data1/bfernando/DE_project/results/qc/Col-0_barcode_rank_forceN.png) has **no knee** —
   a smooth diagonal from rank 1 to ~10^5 with no cliff separating cells from background.
   Per 10x, the cell estimate is unreliable in that case, so there is no data-driven cell
   count. Tested imposing one with `--force-cells` (2000 and 3000) to measure what it buys:

| run                  | cells | med UMI | med genes | med % organellar | total genes |
| -------------------- | ----- | ------- | --------- | ---------------- | ----------- |
| default (EmptyDrops) | 864   | 1,379   | 997       | 3.8              | 20,547      |
| force-cells 2000     | 2,000 | 1,034   | 783       | 5.3              | 20,939      |
| force-cells 3000     | 3,000 | 772     | 611       | 5.3              | 21,095      |

Commands used (`reanalyze`, not `count` — reuses the raw matrix, so no re-alignment):

```bash
cellranger reanalyze --id=Col-0_force2000 \
  --matrix=/data1/bfernando/DE_project/Col-0_v4/outs/raw_feature_bc_matrix.h5 \
  --force-cells=2000 --localcores=8 --localmem=32
```

**CONCLUSION: `--force-cells` was REJECTED — kept the default 864 cells.**
Tripling the cell count (864 -> 3000) added only **548 genes (+2.7%)** while median UMI fell
**44%** (1,379 -> 772), median genes fell **39%** (997 -> 611), and organellar rose 3.8% -> 5.3%.
The extra barcodes contain the same ambient RNA the called cells already sample — more rows in
the matrix, not more information. `--force-cells N` just takes the top N barcodes by total UMI;
it imposes a number rather than discovering cells.

Root cause is upstream, not a parameter: **23% reads in cells** (10x flags <70% as high ambient)
and no knee point at low intact-nuclei recovery / high ambient RNA — a library-prep issue.

**CELL RANGER SHOWED HIGH MULTI MAPPING, SO RAN WITH ALVIN FRY.** (/data1/bfernando/DE_project/notes/01_fromcellrangeroutput.md)

**AlvinFry Runs:**

env at : `/data1/bfernando/conda_envs/alevinfry`

**Recoverable** - reads discarded over  *gene-level ambiguity* . alevin-fry's `cr-like-em` resolves these by EM (** Expectation-Maximization.**) instead of discarding them.

**Not recoverable** - reads that are copies of rRNA. There is no gene-expression information in the molecule. No quantifier retrieves it.

1. **Build the splici reference using `pyroe`**
   alevin-fry matches reads directly against transcript sequences (not against the genome),
   so the transcript sequences have to be extracted first. Cell Ranger's `mkref` output
   cannot be reused for this — it indexes the genome, not the transcriptome.
   1: `make-splici` --> makes the reference transcript sequences. (**splici reference — Spliced + Intronic transcriptome.** Note: this is the *reference* (the sequences). The *index* is built from it separately by `piscem build`.)
   [ because this is single-*nucleus* data. Nuclear RNA is largely unspliced pre-mRNA still carrying introns. A mature-transcript-only reference would discard most of it. Splici = **spli**ced + **i**ntronic — both forms present.]

   **PREPROCESSING OF GTF FILE :**

   1. Cleaned the GTF because `pyroe` (pyranges) crashed on it: `'float' object has no attribute 'rstrip'`. Cause: **5 lines contain no tab characters** — 4 AGAT header lines
      plus 1 blank line — sitting ~640,000 lines deep in the file, because the merge was done
      with `cat araport.gtf contig.gtf`, which put the contig file's header block in the
      middle instead of the top. Those lines parse as `NaN` where a string is expected.
      `cellranger mkref` tolerated them; pyranges does not.
      The awk also renames AGAT's `mRNA` to gffread's `transcript` (the two halves of the
      merged GTF use different vocabularies) and keeps only gene/transcript/exon.
   2. ```
      awk 'BEGIN{FS=OFS="\t"}
           /^#/ {next}
           NF<9 {next}
           $3=="mRNA" {$3="transcript"}
           $3=="gene" || $3=="transcript" || $3=="exon" {print}' \
        reference/build/Athaliana.Col-0.HPIv02.plus_pFusionMYB41.gtf \
        > alvinfry/splici_input.gtf
      ```

**MAKE SPLICI REFERENCE :**

```
pyroe make-splici \
  reference/build/Athaliana.Col-0.HPIv02.plus_pFusionMYB41.fasta \
  alvinfry/splici_input.gtf \
  90 \
alvinfry/splici_HPIv02_pFusionMYB41
```

`90` = R2 length (used to size the flanks added around intron boundaries, so reads straddling
an exon-intron junction still have somewhere to map). Default flank trim is 5, hence `fl85`.

**OUTPUT FILE :**

**`splici_fl85.fa`** --> 98,998 sequences (48,458 spliced + 50,540 intronic). Feeds `piscem build`

**`splici_fl85_t2g_3col.tsv`** — maps each sequence to its gene and splice status. Feeds `alevin-fry quant -m`, for counts table

**PISCEM INDEX BUILDING :**

(Bioconda ships salmon 2.4.1, which removed single-cell mapping. so used piscem based on the alvin fry website )

( its the **indexing** step , it turns your splici FASTA into something  reads can be searched against )

It chops every sequence into **k-mers** — overlapping 31 bp windows — and assembles them into a
**compacted de Bruijn graph**, then puts a hash index (`sshash`) over that graph.

The compaction matters because the 98,998 sequences overlap heavily: isoforms of one gene share
exons, and each gene's intronic record sits adjacent to its spliced one. Rather than storing all
that redundantly, the graph stores each shared stretch **once** as a **unitig** (a non-branching
path) and records which transcripts run through it. That is how ~127 MB of FASTA becomes
85.3 Mb of unitig sequence (620,779 unitigs, avg 137 bp, longest 10,013 bp, k=31).

Then when `map-sc` takes a read, it looks up the read's 31-mers, finds which unitigs they hit,
and reads off which transcripts pass through those unitigs.

Built once and reused — the index depends only on the reference, not the reads, so the same
index serves Col-0, pFACT and pHORST. Index formats are tool-specific: Cell Ranger's
`mkref` output cannot be substituted here.

```
piscem build \
  -s alvinfry/splici_HPIv02_pFusionMYB41/splici_fl85.fa \   # input sequences
  -o alvinfry/piscem_splici_idx \                            # output STEM
  -t 16 \                                                    # threads
  -w alvinfry/piscem_workdir                                 # scratch space
```

**OUTPUT :**

```
piscem_splici_idx.ssi        65 MB   sshash k-mer index
piscem_splici_idx.ssi.mphf    3 MB   minimal perfect hash
piscem_splici_idx.tdct       1.6 GB  dictionary (the bulk)
piscem_splici_idx.tct        17 MB   contig table
piscem_splici_idx.ctab        6 MB   contig -> reference mapping
piscem_splici_idx.ectab       3 MB   equivalence classes
piscem_splici_idx.refinfo     4 MB   reference names and lengths
```

**MAP READS TO RAD :**

```
piscem map-sc \
  -i alvinfry/piscem_splici_idx \
  -g '1{b[16]u[12]x:}2{r:}' \
  -1 /data2/nhartwick/law_lab_myb41_scrna/rawreads/Col-0_S1_L003_R1_001.fastq.gz \
  -2 /data2/nhartwick/law_lab_myb41_scrna/rawreads/Col-0_S1_L003_R2_001.fastq.gz \
  -o alvinfry/Col-0_map -t 16
```

```
Geometry 1{b[16]u[12]x:}2{r:} = from R1 take 16 bp barcode + 12 bp UMI, discard the rest; R2 is the cDNA. This is the 10x 3′ v4 layout.
```

Confirmed from Piscem logs :

```
Protocol: custom (bc_len=16, umi_len=12)
```

**RAD** = Reduced Alignment Data. Not a BAM — no sequences, no quality scores, no coordinates. Only the mapping outcome per read, which is why it's compact and everything downstream runs in minutes.

**What this step does:** for every read pair — pull barcode + UMI off R1, match R2 against the index, write one RAD record (barcode, UMI, list of sequences hit).

**What it deliberately does not do:** no counting, no UMI deduplication, no barcode validation, **and no multimapping resolution** — if a read hits five sequences, all five are recorded. Cell Ranger discards multimappers at alignment time and they are gone forever; here the ambiguity survives into the RAD so the resolution strategy can be chosen later.

OUTPUT :

```
Mapped 88187028/458636011 reads (19.23%), 0 poisoned, 243.5s
```

Output: `alvinfry/Col-0_map/map.rad` (1.68 GB) + `map_info.json`

---

**TWO PIPELINES ARE CONSISTENT**

From alevin-fry, only **19.23%** of reads match the transcriptome.

So **80.8% of the library is not transcriptomic at all.**

This is an independent measurement of the same thing the Cell Ranger
metrics implied, arrived at by completely different means.

|                  | Cell Ranger                                              | alevin-fry / piscem                                 |
| ---------------- | -------------------------------------------------------- | --------------------------------------------------- |
| total reads      | 458,636,011                                              | 458,636,011                                         |
| what it measures | reads aligning**confidently to the transcriptome** | reads**matching splici transcript sequences** |
| result           | **16.2%**                                          | **19.23%**                                    |

Same library (read totals identical), two different aligners, two different reference types (genome + GTF vs. splici transcriptome), two different algorithms  and both say **only ~16-19% of this library is usable transcriptomic signal.**

The ~3 percentage-point difference is the multimapping recovery that piscem retains reads that map to several transcripts, which Cell Ranger discards at alignment time.

**WHAT THE SWITCH TO ALEVIN-FRY ACTUALLY GAINED:**

|                                            | usable reads   |
| ------------------------------------------ | -------------- |
| Cell Ranger (confidently -> transcriptome) | 16.2% = 74.3M  |
| alevin-fry + splici                        | 19.23% = 88.2M |

**+13.9M reads, ~19% more usable data.**

**CHECK FOR THE 80% OF MISSING READS :**

env : `conda activate de_project`

```bash
samtools idxstats /data1/bfernando/DE_project/Col-0_v4/outs/possorted_genome_bam.bam
```

Columns are: reference, length (bp), mapped reads, unmapped reads.

```
Chr1            30427671    18338962      0
Chr2            19698289   146114338      0
Chr3            23459830   217895296      0
Chr4            18585056    10835887      0
Chr5            26975502    16438079      0
ChrC              154478     3496063      0
ChrM              366924    36395838      0
pFusionMYB41        1866          67      0
*                      0           0    9121481
```

**WITH READS/KB per CHR**

```
samtools idxstats /data1/bfernando/DE_project/Col-0_v4/outs/possorted_genome_bam.bam \
  | awk 'BEGIN{OFS="\t"; print "ref","len_kb","reads","reads_per_kb"}
         $1!="*"{printf "%s\t%.1f\t%d\t%.1f\n", $1, $2/1000, $3, $3/($2/1000)}
         $1=="*"{printf "unmapped\t-\t%d\t-\n", $4}'
```

```
ref     len_kb  reads   reads_per_kb
Chr1    30427.7 18338962        602.7
Chr2    19698.3 146114338       7417.6
Chr3    23459.8 217895296       9288.0
Chr4    18585.1 10835887        583.0
Chr5    26975.5 16438079        609.4
ChrC    154.5   3496063         22631.5
ChrM    366.9   36395838        99191.8
pFusionMYB41    1.9     67      35.9
unmapped        -       9121481 -
```

**METHOD — why normalise by length**

Raw read counts mean nothing on their own; a longer chromosome should collect more reads.
Dividing by length is what makes the anomaly visible.

**Chr1, Chr4 and Chr5 agree at ~600 reads/kb** (602.7, 583.0, 609.4 — within 5% of each other).
Three independent chromosomes converging is the empirical baseline for "normal". Everything
else is measured against it. Baseline used = **600.3 reads/kb** (pooled Chr1+Chr4+Chr5).

| ref            | reads/kb | vs baseline      | % of mapped reads |
| -------------- | -------- | ---------------- | ----------------- |
| Chr1           | 602.7    | 1.0x             | 4.1%              |
| **Chr2** | 7,417.6  | **12.4x**  | **32.5%**   |
| **Chr3** | 9,288.0  | **15.5x**  | **48.5%**   |
| Chr4           | 583.0    | 1.0x             | 2.4%              |
| Chr5           | 609.4    | 1.0x             | 3.7%              |
| **ChrC** | 22,631.5 | **37.8x**  | 0.8%              |
| **ChrM** | 99,191.8 | **165.8x** | **8.1%**    |
| pFusionMYB41   | 35.9     | 0.06x            | ~0%               |

**Chr2 + Chr3 + ChrM alone hold 89% of every mapped read.**

**EXPECTED at baseline (length_kb x 600.3):**

```
Chr2 : 19698.3 x 600.3 = 11.8M
Chr3 : 23459.8 x 600.3 = 14.1M
ChrM :   366.9 x 600.3 =  0.22M
ChrC :   154.5 x 600.3 =  0.09M
```

**EXCESS over baseline:**

```
Chr2   146.1M observed  -  11.8M expected  =  134.3M excess
Chr3   217.9M observed  -  14.1M expected  =  203.8M excess
ChrC     3.5M observed  -   0.09M expected =    3.4M excess
ChrM    36.4M observed  -   0.22M expected =   36.2M excess
                                              ─────────────
                                               377.7M excess

377.7M / 458,636,011  =  82.3% of the library
```

**CONFIRMATION — two independent methods agree:**

| method                                    | evidence                             | "not usable"    |
| ----------------------------------------- | ------------------------------------ | --------------- |
| `samtools idxstats` density subtraction | where reads land on chromosomes      | **82.3%** |
| piscem splici mapping                     | do reads match annotated transcripts | **80.8%** |

Within 1.5 percentage points, from completely different evidence, on the same library
(read totals identical: 458,636,011).
