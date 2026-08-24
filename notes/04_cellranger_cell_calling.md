# How Cell Ranger counts cells (cell-calling) — notes

**Purpose:** reference for how `cellranger count` turns reads into "cells," why our Col-0
run recovered few cells, and what `--force-cells` does. Verified against the 10x docs
(links at bottom) and our own Col-0_v4 numbers.

---

## 1. The hierarchy: reads → barcodes → UMIs → genes

Every read carries three pieces of information:

| From          | Piece                        | Identifies                                                         |
| ------------- | ---------------------------- | ------------------------------------------------------------------ |
| R1, bp 1–16  | **cell barcode**       | which droplet / nucleus                                            |
| R1, bp 17–28 | **UMI**                | which individual RNA molecule                                      |
| R2 (cDNA)     | **alignment position** | which**gene** (R2 maps to the genome → falls in a GTF gene) |

- **UMI count for a barcode** = number of *distinct* molecules captured (PCR duplicates —
  same barcode+UMI+gene — are collapsed to 1). It is a **count of molecules, not a score.**
- **Real cell** = intact nucleus = lots of RNA = **high** UMI count.
- **Empty droplet (background)** = no nucleus, only a little ambient RNA = **low** UMI count.

Two output matrices (genes × barcodes):

- `raw_feature_bc_matrix/` = **all** barcodes (cells + background). Col-0: 434,353 barcodes.
- `filtered_feature_bc_matrix/` = **only called cells**. Col-0: 864 cells.
- The difference between them **is** the cell-calling result.

> Counting distinct R1 barcodes does NOT give the cell count — it gives millions
> (empty droplets + sequencing errors). Cell-calling is what separates real cells out.

---

## 2. The barcode-rank plot (the key QC picture)

Sort barcodes by UMI count (high→low) and plot on log-log axes. **Each y-axis gridline =
one order of magnitude (10×).**

- **Good sample:** a tight high band of cells (spanning ~1 order of magnitude), then a
  **steep cliff** dropping to a much lower background band. The cliff = clean separation.
- **Bad sample (ours):** a **smooth gradual slope**, no cliff — cells and background smear
  together across several decades. Per 10x: *"the estimated number of cells is unreliable
  because the algorithm has trouble discerning cells from background."*

---

## 3. The cell-calling algorithm (default, automatic)

Two data-driven steps using both **quantity (UMIs)** and **quality (gene profile)**:

### Step 1 — OrdMag (Order of Magnitude)

Automatically picks a UMI cutoff, assuming real cells cluster within ~1 order of magnitude:

1. Take the top N barcodes by UMI (N = expected cells).
2. `m` = **99th-percentile** UMI among those top N (a robustly-high "typical strong cell"
   value; 99th percentile, not the max, to ignore doublet/outlier spikes).
3. Cutoff = **`m / 10`** (one order of magnitude below).
4. Call every barcode with UMI > m/10 a cell.

Works only when cells form a ~10× band sitting clearly above a background band with a gap.

**Call every barcode with UMI count > m/10 a cell, where m = the 99th-percentile UMI value.**

**OrdMag = "cell if UMI > m/10, with m = the 99th-percentile UMI of the top barcodes."**

### Step 2 — EmptyDrops

Refines using the **gene profile**, not just UMI count:

1. Build an **ambient/background model** (a multinomial over the ambient gene profile,
   estimated from low-UMI barcodes via Simple Good-Turing smoothing).
2. For **each candidate barcode (>500 UMIs, above ambient level, not already an OrdMag cell),**
   test: **does its gene profile differ significantly from the ambient model?**
   - **Yes** → real cell (even at modest UMI count) → **rescued**.
   - **No** → ambient soup (even at high UMI count) → **rejected**.

**Consequence:** the called-cell set is **NOT simply "top-N by UMI."** EmptyDrops can
reject high-UMI ambient barcodes and rescue lower-UMI barcodes with distinct profiles.

---

## 4. `--force-cells N` vs automatic cell-calling

|                            | Auto (OrdMag + EmptyDrops)  | `--force-cells N`                |
| -------------------------- | --------------------------- | ---------------------------------- |
| Who sets the number        | the data (algorithm)        | **you**                      |
| Signals used               | UMI count**+ gene profile** | **UMI count only**           |
| Reject high-UMI ambient?   | ✅ yes                      | ❌ no (if in top N, it's a "cell") |
| Rescue low-UMI real cells? | ✅ yes                      | ❌ no (hard cutoff at rank N)      |

`--force-cells N` = *"ignore the algorithm; rank barcodes by UMI; keep the top N."* It
**imposes** a number, it does not **discover** cells. Fast way to apply it (no re-align):
`cellranger reanalyze --matrix raw_feature_bc_matrix.h5 --force-cells N`.
Use only with outside knowledge (loaded-cell count) or a hopeless knee — and always vet
the added barcodes' quality (they may be ambient).

Tunable-parameter note: `cellranger count` exposes only `--expect-cells` (soft hint to
OrdMag) and `--force-cells` (hard override). Internal knobs (m/10, EmptyDrops FDR/lower
bound) are NOT CLI-tunable — for that, run `DropletUtils::emptyDrops` (R) on the raw matrix.

---

## 5. Key QC metrics (and 10x thresholds)

- **Estimated Number of Cells** — the cell-calling result.
- **Fraction Reads in Cells** — fraction of usable reads in cell barcodes. Per 10x,
  **< 70% indicates high ambient RNA** partitioned into all GEMs.
- **Barcode-rank knee** — a steep cliff = good separation; no cliff = unreliable cell count.
- High **mitochondrial/organellar** fraction → dead/dying cells / ambient.

---

## 6. Col-0_v4 case study (why we recovered few cells)

- **864 cells**; median 1,380 UMI/cell, ~1,000 genes/cell (reasonable *for nuclei*).
- **Fraction reads in cells = 23%** (< 70% → high ambient, per 10x).
- **No barcode-rank knee** — smooth slope.
- **Overlap (no cliff), measured:** lowest called cell = **502 UMIs**, but the highest
  **background** (non-cell) barcode = **5,153 UMIs** — an empty droplet out-counts most
  cells. UMI count alone cannot separate cells from background here.
- Root cause: **low intact-nuclei recovery / high ambient RNA** (library issue), plus
  ~20–25% TSO artifact reads. NOT a reference or mapping problem (reference validated;
  98% reads map to genome; 20,547 genes detected; `AT4G28110.Fusion` = 0 in WT control).
- `--force-cells` can raise the number, but the added barcodes have steadily lower
  UMI/genes (865–2000: ~799 UMI; 3001–5000: ~362 UMI) — top-by-UMI barcodes EmptyDrops
  did **not** endorse. Organellar % stays ~5% (not obvious dead-cell junk, but ambiguous).

---

## Sources (10x Genomics)

- Gene Expression Algorithm (OrdMag + EmptyDrops): https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm
- Barcode Rank Plot: https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-barcode-rank-plot
- Quality Assessment (Fraction Reads in Cells < 70% → ambient): https://www.10xgenomics.com/analysis-guides/quality-assessment-using-the-cell-ranger-web-summary
- Introduction to Ambient RNA Correction: https://www.10xgenomics.com/analysis-guides/introduction-to-ambient-rna-correction
