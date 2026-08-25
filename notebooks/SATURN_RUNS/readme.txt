# SATURN: integrating 4 Arabidopsis datasets

SATURN does not specify if we need to use the raw counts or normalized counts, 
in our runs we have used COUNT values in .X, based on reasoning that pretraining uses an scVI ZINB loss
and HVG selection runs flavor='seurat_v3', which assumes counts.

Prep is in scripts/SATURN_RUN_BASED_SCRIPTS/dataprep.py
Converting each data to its own un-normalozed raw counts form.

Paths below are relative to two configurable roots:
    $DE_PROJECT_DIR   this repository
    $PENNYCRESS_DIR   directory holding adata_arabidopsis_{shahan,nucleus}/


## 1) Shahan            110,427 x 28,958

    Already raw counts -- dtype int64.
    NO CHANGE NEEDED.

    $PENNYCRESS_DIR/adata_arabidopsis_shahan/shahan_full_devlabel.h5ad


## 2) Arabidopsis nucleus   19,755 x 31,559

    NOT counts. float64, values 0.7326..7.0149, a max in single digits means
    log space, since log1p compresses 5124 down to 8.5.

    Identified the exact transform by inverting it: after np.expm1() every cell
    summed to exactly 10000.00 (CV = 0.000000), i.e.
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    Counts recovered without the source file: expm1 to undo the log, then divide
    each row by its own smallest nonzero value. 

    Saved to:
    $PENNYCRESS_DIR/adata_arabidopsis_nucleus/arabidopsis_ch_snRNA_devlabel_RAWCOUNTS_.h5ad

    Also note: 31,559 genes is larger than Arabidopsis's ~27,655 protein-coding
    genes because this annotation includes non-coding loci. 
    3,982 of them have no ESM3 embedding and SATURN drops them at load, leaving 27,577.


## 3) MYB41 transgene samples        in $DE_PROJECT_DIR/alvinfry/h5ad/

    Counts, but FRACTIONAL -- alevin-fry cr-like-em splits an ambiguous UMI probabilistically across genes, 
    so values like 0.01 appear. 
    Not normalization:
        per-cell totals vary 500..157,206 (CV ~2.0-2.5), and most values are already
        whole. Fixed by np.round() + eliminate_zeros() + astype(int32).

    Col-0   -> Col-0_INT_COUNTS.h5ad    (2442, 27656)
    pFACT   -> pFACT_INTCOUNTS.h5ad     (3304, 27656)
    pHORST  -> pHORST_COUNTS.h5ad       (2937, 27656)

    Transgene check after rounding -- cell counts must hold at 0 / 198 / 34.
    Col-0's total drops to 0 because its 2.00 was entirely sub-1 fractional dribs,
    which is why the detection threshold is >= 1 UMI rather than > 0.

