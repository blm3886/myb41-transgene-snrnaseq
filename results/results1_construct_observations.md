# Results 1 — Observations on the MYB41 transgene constructs & custom reference

**Project:** snRNA-seq transgene detection in *Arabidopsis* (Col-0)
**Constructs:** `pDEST_pFACT_MYB41_3UTR` and `pDEST_pHORST_MYB41_3UTR`
**Source notebooks:** [`notebooks/01_compare_constructs_2.ipynb`](../notebooks/01_compare_constructs_2.ipynb),
[`notebooks/02_build_reference.ipynb`](../notebooks/02_build_reference.ipynb)

---

## Summary

We have compared the two transgenic constructs directly from their plasmid maps and built a
single custom reference to detect transgene expression. Based on this:

1. **FACT and HORST are identical except for the promoter** : I have done base-by-base comparison from the raw DNA.
2. **A transcribed cloning scar ("Scar 2") is the detection marker** : As described in the meeting ,This is too verified to be identical in both constructs and absent from the native genome, so reads over it uniquely mark the transgene.
3. **Thus, One reference suffices for both constructs**; FACT vs HORST lines are distinguished by
   sample metadata, not by reads.

---

## 1. The two constructs differ *only* in the promoter

Both plasmids share the same cassette layout (T-DNA between the Left/Rright borders). Comparing the full sequences:

|                                | FACT                                              | HORST                                                  |
| ------------------------------ | ------------------------------------------------- | ------------------------------------------------------ |
| Plasmid size                   | 11,506 bp                                         | 12,006 bp                                              |
| Identical prefix (from pos 1)  | 9,293 bp                                          | 9,293 bp                                               |
| Identical suffix (to the end)  | 628 bp                                            | 628 bp                                                 |
| **Sole differing block** | pFACT promoter,**1,582 bp** (9,295–10,876) | proCYP86A1 promoter,**2,082 bp** (9,295–11,376) |

- The two constructs are **byte-for-byte identical** across 9,293 bp of prefix + 628 bp of
  suffix. They diverge in **exactly one contiguous block: the promoter.**
- The **500 bp** length difference between the two constructs:
  - The promoter is also the only region that changes  **length**
  - HORST's promoter is 500 bp longer than FACT's (2,082 vs 1,582 bp)
  - This account for the 500 bp difference in total plasmid size (12,006 vs 11,506 bp)

---

## 2. Scar 2 is the transgene detection marker

The cassette contains two cloning scars. On the minus strand (transcription runs high→low coordinate):

| Scar             | Plasmid coords | Length | Location                          | Transcribed?  | In snRNA-seq?               |
| ---------------- | -------------- | ------ | --------------------------------- | ------------- | --------------------------- |
| **Scar 2** | 8,024–8,051   | 28 bp  | between 3′UTR and MYB41 gene end | **Yes** | **Yes — the marker** |
| Scar 1           | 9,268–9,294   | 27 bp  | between MYB41 gene and promoter   | No            | No                          |

**Scar 2 sequence:** `AGGCCACTTTGTACAAGAAAGCTGGGTC`

Three independent checks confirm Scar 2 is a valid, specific marker:

- **Identical in both constructs** — the 28 bp sequence is byte-for-byte the same in FACT and
  HORST (Scar 2 falls inside the shared, identical region).
- **Absent from the native genome** — 0 occurrences (both strands) of the 28 bp across the
  entire Col-0 HPIv02 genome (7 sequences). Reads matching Scar 2 can therefore *only*
  originate from the transgene.
- **Correctly placed in the reference** — Scar 2 maps to contig positions **390–417**,
  consistent with the annotation's transcribed `insertion` feature (388–417).

*(Scar 1 differs by 1 bp between constructs at its promoter-facing end, but it is untranscribed
and irrelevant to detection.)*

---

## 3. Additional observations

- **MYB41 = AT4G28110.** The native gene sits on **Chr4, minus strand, 13,968,029–13,969,525**
  (Araport11). The transgene copy is annotated separately as **`AT4G28110.Fusion`**, so native
  and transgene expression are distinguishable.
- **Both transgenes are on the minus strand** (promoter, scars, MYB41, UTR all `−`).
- **Verified About the Provided contig:** the `pFusionMYB41` contig (1,866 bp) is an **exact copy of the FACT plasmid, positions 7,635–9,500** (found only in FACT — it carries a ~206 bp pFACT
  promoter tail — but its transcribed region is identical in HORST, so it detects either). : [`01_compare_constructs_2.ipynb`](vscode-webview://04dbh5n34pr2pfrdojmkloou61mvj8kerhcrq5mhba994k7sc0l8/notebooks/01_compare_constructs_2.ipynb)
- The contig annotation (GTF) carries **both  isoform : `AT4G28110.Fusion.1` **and** `AT4G28110.Fusion.2`. : [`02_build_reference.ipynb`](vscode-webview://04dbh5n34pr2pfrdojmkloou61mvj8kerhcrq5mhba994k7sc0l8/notebooks/02_build_reference.ipynb)**

---

## 4. The custom reference (built in notebook 02)

A single reference was assembled for aligning all samples:

- **Base genome:** Col-0 **HPIv02** (Salk HPI). *(HPIv01 and HPIv02 have byte-identical genome and gene models; v02 chosen.)*
- **Combined genome FASTA:** Col-0 chromosomes **+ the transgene contig = 8 sequences**
  (Chr1–5, ChrC, ChrM, `pFusionMYB41`). The contig is added as a standalone sequence — the transgene is not inserted into the chromosome.
- **Combined annotation:** base GFF3 converted to GTF and merged with the contig GTF; both the native and `AT4G28110.Fusion` gene models are present.
- **Validated:** transgene contig length 1,866 bp; Scar 2 at contig 390–417 matches the expected sequence; native + fusion MYB41 both annotated; `samtools faidx` succeeds.

### Build notes (fixes)

- The provided contig FASTA was **malformed** (no `>` header) and was rewritten to valid FASTA.
- The source GFF3 has **no `exon` features** (only gene/mRNA/CDS/UTR); `gffread -T` left ~3,824
  CDS-only transcripts without exons, which Cell Ranger `mkref` rejects. Fixed by re-converting
  with **`gffread --force-exons`** (0 transcripts without exons).

---

## 5. Conclusion & implications for the analysis

- **Build one reference, not two.** The constructs share an identical detectable sequence, so a single spike-in works for both. **FACT vs HORST lines are told apart by sample metadata**
- **Detect the transgene via Scar 2.** Reads mapping uniquely to the `pFusionMYB41` contig (Scar 2, contig 390–417) mark transgene expression. These should appear **only in the
  pFACT/pHORST samples, not in the Col-0 wild-type control**, and can be projected onto the cell clusters to localize which nuclei express the transgene.

---

*Reproducibility: all numbers above are computed and asserted in the two source notebooks,
which run against the plasmid maps and the Col-0 HPIv02 reference.*
