
# Cell calling on the alevin-fry matrices

---

The counts matrices are **barcode x gene**, not **cell x gene**.

```
barcode  ->  gel bead  ->  droplet (GEM)  ->  0, 1, or 2+ nuclei
```

A barcode identifies a **droplet**, not a cell. Most droplets never caught a nucleus, but they still get a barcode, and they still capture free-floating "ambient" RNA from nuclei that broke during prep. So an empty droplet produces reads, gets counted, and appears as a row in the matrix.

Nothing in the pipeline so far has tried to separate them:

- `map-sc` recorded whatever 16 bases sat at the start of R1., barcode or droplet.
- `generate-permit-list --unfiltered-pl` asked *"*is this a real 10x barcode?" — a manufacturing
  question, answered against the whitelist. Not *"did this droplet contain a nucleus?"*

`--unfiltered-pl` **deliberately** skips cell calling: it keeps every barcode with
`--min-reads` >= 10 (default) and leaves the decision to downstream analysis.

**Current barcode counts:**

| sample | barcodes in matrix |
| ------ | ------------------ |
| Col-0  | 192,259            |
| pFACT  | 200,283            |
| pHORST | 193,477            |



For comparison, Cell Ranger on the same Col-0 library:

```
Col-0_v4/outs/raw_feature_bc_matrix/        434,353 barcodes   <- comparable to ours
Col-0_v4/outs/filtered_feature_bc_matrix/       864 cells      <- after cell calling
```

**The difference between those two Cell Ranger files IS the cell-calling result.**
alevin-fry only produces the equivalent of the first one.

---

## 2. ATTEMPT 1 — KNEE METHOD (`generate-permit-list -k`)

`generate-permit-list` has five mutually exclusive modes:

```
<--knee-distance | --expect-cells | --force-cells | --valid-bc | --unfiltered-pl>
```

We originally used `-u/--unfiltered-pl`. The alevin-fry tutorial
(https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/)
uses `-k/--knee-distance`, which *does* call cells. Tested it on Col-0:

```bash
alevin-fry generate-permit-list -i alvinfry/Col-0_map -d fw -k \
  -o alvinfry/Col-0_quant_knee
alevin-fry collate -i alvinfry/Col-0_quant_knee -r alvinfry/Col-0_map -t 16
alevin-fry quant -i alvinfry/Col-0_quant_knee \
  -m alvinfry/splici_HPIv02_pFusionMYB41/splici_fl85_t2g_3col.tsv \
  -t 16 -r cr-like -o alvinfry/Col-0_quant_knee
```

`generate_permit_list.json` confirms `"fmeth": "KneeFinding"`, `"permit-list-type": "filtered"`.

### RESULT — the knee method FAILEDThe two distributions are almost identical:

```
                      barcodes    totalUMI    median  mean   max
knee (-k)              152,757   17,738,216      68    116   65,273
unfiltered (-u)        192,259   17,592,435      62     92   65,337

                      >100 UMI   >500 UMI   >1000 UMI
knee (-k)               29,311      2,961       1,174
unfiltered (-u)         28,568      2,934       1,162
```

The knee filter dropped ~39,500 low-count barcodes and left everything else untouched.

### WHY IT FAILED

1. **There is no knee.** The barcode-rank curve is a smooth diagonal from rank 1 to ~10^5 with
   no cliff (see `results/qc/Col-0_barcode_rank_forceN.png`). A knee-finding algorithm run on a
   curve with no knee returns a number, but it is not reading a real feature.
2. **`-k` works on read counts, not UMI counts** — a cruder signal than the UMI-based curve,
   applied before deduplication.
3. **Count alone cannot separate cells from background in this library.** Measured overlap from
   the Cell Ranger run: lowest *called* cell = **502 UMIs**, highest *rejected* barcode =
   **5,153 UMIs**. An empty droplet out-counted most real cells tenfold.

**This is a result, not just a failed attempt:** two independent tools (Cell Ranger's
barcode-rank plot and alevin-fry's knee finder) now confirm the same pathology. Cell calling
here must use a **profile-based** method, not a count-threshold method.

> `alvinfry/Col-0_quant_knee/` is kept as evidence but **must not be used downstream.**
> The runs of record remain `Col-0_quant`, `pFACT-MYB41_quant`, `pHORST-MYB41_quant`.

---

## 3. WHERE THE REAL CELLS PROBABLY ARE

Even though thresholding cannot *call* cells here, it brackets the plausible answer:

| threshold  | barcodes (Col-0) |
| ---------- | ---------------- |
| > 100 UMI  | 28,568           |
| > 500 UMI  | **2,934**  |
| > 1000 UMI | **1,162**  |

Cell Ranger called **864** cells with a minimum of 502 UMI, which sits between these. So the
true count is plausibly **~800-1,500**, not 152,757. Use this as a sanity check on whatever
emptyDrops returns.

---

## 4. ATTEMPT 2 — emptyDrops (PLANNED)

Thresholding on UMI count is exactly what fails here. **emptyDrops tests gene composition
instead of quantity:**

1. Take every barcode with UMI <= `lower` (default 100); assume these are definitely empty.
2. Pool their counts into an **ambient profile** — a probability distribution over genes.
3. For each candidate barcode, ask: given its total UMI count, how surprising is its observed
   gene composition under that ambient distribution? (multinomial, p-value by Monte Carlo)
4. Low p -> composition too distinctive to be soup -> **real cell**. FDR-corrected.
5. Very high-UMI barcodes are auto-retained regardless.

This is why it can keep a 502-UMI barcode with a distinctive profile while rejecting a
5,153-UMI barcode that just looks like a large scoop of average soup. It is also the same
algorithm Cell Ranger ran internally to get 864 — so Col-0 gives us a validation target.

```bash
conda install -n r-env -c conda-forge -c bioconda bioconductor-dropletutils
```

Export from Python (matrix must be U+S+A summed, i.e. 27,656 genes, not 82,968 columns):

```python
import pyroe, numpy as np
from scipy.io import mmwrite
SAMPLES = {"Col-0":"Col-0", "pFACT":"pFACT-MYB41", "pHORST":"pHORST-MYB41"}
for name, d in SAMPLES.items():
    a = pyroe.load_fry(f"alvinfry/{d}_quant", output_format="snRNA")
    mmwrite(f"alvinfry/celcall_{name}.mtx", a.X.T)              # genes x barcodes
    np.savetxt(f"alvinfry/celcall_{name}_bc.txt", a.obs_names.values, fmt="%s")
```

```r
library(DropletUtils); library(Matrix)
for (s in c("Col-0","pFACT","pHORST")) {
  m  <- readMM(sprintf("alvinfry/celcall_%s.mtx", s))
  bc <- readLines(sprintf("alvinfry/celcall_%s_bc.txt", s))
  set.seed(42)
  e    <- emptyDrops(m, lower = 100)
  keep <- which(!is.na(e$FDR) & e$FDR <= 0.01)
  cat(s, "cells called:", length(keep), "\n")
  writeLines(bc[keep], sprintf("alvinfry/cells_%s.txt", s))
}
```

Knobs: **`lower`** (100 default) sets which barcodes build the ambient profile — try 100 and 200
and see how much the answer moves. **FDR** 0.01 is standard, 0.001 stricter.

**VALIDATION:** Col-0 must land near 864 (anything ~800-1,500 is credible). If it returns tens
of thousands, the ambient profile estimate is the problem and CellBender becomes the next step
rather than a refinement.

---

## 5. CELL CALLING IS NOT AMBIENT REMOVAL

Two distinct problems, and this library has both:

|                           | question                                             | granularity             | tool               |
| ------------------------- | ---------------------------------------------------- | ----------------------- | ------------------ |
| **Cell calling**    | is this droplet a cell?                              | per**barcode**    | emptyDrops         |
| **Ambient removal** | of the counts in this cell, how many came from soup? | per**gene, cell** | CellBender / SoupX |

emptyDrops decides which barcodes to keep — it does not clean the ones it keeps. With
**23% reads in cells** (10x flags <70% as high ambient), the retained cells will still carry
substantial contamination.

### Why this matters for the transgene specifically

**In pFACT, the ambient soup contains transgene transcript.** Consequences both ways:

- *Helpful:* emptyDrops builds its ambient profile including the transgene, so a droplet is not
  promoted to "cell" merely for carrying transgene reads.
- *Harmful:* every real cell in pFACT picks up some transgene counts from soup, whether or not
  that nucleus expressed it. A cell with 1 Fusion UMI is **not** evidence of expression in that
  cell.

**Detection is already settled** and does not depend on any of this — pFACT 1,312 UMIs vs a
wild-type floor of 1, with Col-0 controlling for ambient by construction. But **localisation**
("which cell types express the transgene") cannot be answered from raw per-cell counts, because
soup would paint low-level signal across every cell type.

Expect the 1,182 transgene-positive **barcodes** in pFACT to fall substantially once real cells
are called — most are empty droplets that caught ambient transgene RNA. After cell calling, look
at the *distribution* of per-cell Fusion counts rather than just how many are non-zero; a real
expressing population should sit clearly above background, and Col-0 gives the shape of
ambient-only for comparison.

---

## TO DO

- [ ] Install `bioconductor-dropletutils` in `r-env`
- [ ] Run emptyDrops on all three unfiltered matrices
- [ ] Validate Col-0 against 864 (accept ~800-1,500)
- [ ] Sensitivity check: `lower = 100` vs `200`
- [ ] QC the called cells (UMI, genes, organellar %) — if they look like soup, cell calling
  alone was not enough
- [ ] Re-count transgene-positive **cells** (not barcodes) per sample
- [ ] Decide whether ambient removal (CellBender) is required before localisation
