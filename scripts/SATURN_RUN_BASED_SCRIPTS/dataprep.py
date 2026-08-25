"""Prepare the four Arabidopsis arms for a SATURN run.

SATURN requires raw COUNT values in `.X` (see the snap-stanford/SATURN README:
"scRNA-seq count values in the `X` field"), because pretraining uses an scVI ZINB
loss and highly-variable-gene selection runs `flavor='seurat_v3'`, which assumes
counts. This script brings every arm onto that convention:

  ARM 1  shahan   - already int64 counts, inspected only
  ARM 2  nucleus  - stored as log1p(normalize_total(1e4)); counts are recovered here
  ARM 3  Col-0  ) - alevin-fry `cr-like-em` emits FRACTIONAL counts (an ambiguous UMI
  ARM 4  pFACT  )   is split probabilistically across genes), so these are rounded
         pHORST )   to integers rather than un-normalised

Paths are configurable. Set these before running:

    export DE_PROJECT_DIR=/path/to/DE_project        # optional, inferred from this file
    export PENNYCRESS_DIR=/path/to/charlotte_pennycress
"""

import os
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
from scipy.sparse import issparse

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Configurable roots - override via environment rather than editing this file.
# ---------------------------------------------------------------------------
BASE_DIR = Path(os.environ.get('DE_PROJECT_DIR', Path(__file__).resolve().parents[2]))

_penny = os.environ.get('PENNYCRESS_DIR', '')
if not _penny:
    raise SystemExit(
        'Set PENNYCRESS_DIR to the directory containing '
        'adata_arabidopsis_shahan/ and adata_arabidopsis_nucleus/'
    )
PENNYCRESS_DIR = Path(_penny)

SHAHAN_H5AD  = PENNYCRESS_DIR / 'adata_arabidopsis_shahan'  / 'shahan_full_devlabel.h5ad'
NUCLEUS_H5AD = PENNYCRESS_DIR / 'adata_arabidopsis_nucleus' / 'arabidopsis_ch_snRNA_devlabel.h5ad'
NUCLEUS_OUT  = PENNYCRESS_DIR / 'adata_arabidopsis_nucleus' / 'arabidopsis_ch_snRNA_devlabel_RAWCOUNTS_.h5ad'
MYB          = BASE_DIR / 'alvinfry' / 'h5ad'

# ============================ ARM 1 : SHAHAN ============================
shahan_a = ad.read_h5ad(SHAHAN_H5AD, backed='r')
print(shahan_a.shape) # 110,427 × 28,958
print(shahan_a)


#AnnData object with n_obs × n_vars = 110427 × 28958 backed at '<SHAHAN_H5AD>'
    #obs: 'Celltype', 'celltype_fine', 'time_celltype', 'timezone', 'cell_cycle', 'orig.ident', 'n_counts', 'n_genes', 'dev_label'
    #layers: None (.X)

print(shahan_a.X) # CSRDataset: backend hdf5, shape (110427, 28958), data_dtype int64
print(shahan_a.X._data)
print(shahan_a.X._indices)
print(shahan_a.X._indptr)

# ============================ ARM 2 : NUCLEUS ============================
nuc = ad.read_h5ad(NUCLEUS_H5AD, backed='r')
print(nuc.shape)   # (19755, 31559)
print(nuc)
# AnnData object with n_obs × n_vars = 19755 × 31559 backed at '...'
#     obs: 'Celltype', 'dev_label'
#     layers: None (.X)

print(nuc.X)         # _CSCDataset  <- CSC, not CSR. column-compressed.
print(nuc.X.dtype)   # float64      <- NOT counts: log1p(normalize_total(1e4))
d = nuc.X._data[:]
d   = d[d != 0] 
print(f'{d.min():.4f} / {d.max():.4f}') #### 0.7326 / 7.0149

###################

X = nuc.to_memory().X.tocsr()
X.data = np.expm1(X.data) # Undo the log

for i in range(X.shape[0]):
    lo, hi = X.indptr[i], X.indptr[i+1]
    if hi > lo:    
        X.data[lo:hi] /= X.data[lo:hi].min()  
print(np.allclose(X.data, np.round(X.data), atol=0.01))
result = np.allclose(X.data, np.round(X.data), atol=0.01)

if result:
    # NB: must assign back to X.data. Writing to a new attribute (e.g. X.data_norm)
    # leaves X.data untouched, and the astype below would then TRUNCATE rather than
    # round -- 2.9999998 becomes 2 instead of 3.
    X.data = np.round(X.data)
    X.eliminate_zeros()
    X = X.astype(np.int32)

    nuc_counts = ad.AnnData(X=X, obs=nuc.obs.copy(), var=nuc.var.copy())
    nuc_counts.write_h5ad(NUCLEUS_OUT, compression='gzip')
    print('wrote', NUCLEUS_OUT)
else:
    print('WARNING: count recovery failed the integer check - not writing ARM 2')

print("^^^^^^^^^^^^ARABIDOPSIS NUCLEUS^^^^^^^^^^^^")

# ================== ARMS 3/4 : MYB41 TRANSGENE SAMPLES ==================
col0   = ad.read_h5ad(f'{MYB}/Col-0_cr-like-em_filtered.h5ad')    # (2442, 27656)
pfact  = ad.read_h5ad(f'{MYB}/pFACT_cr-like-em_filtered.h5ad')    # (3304, 27656)
phorst = ad.read_h5ad(f'{MYB}/pHORST_cr-like-em_filtered.h5ad')   # (2937, 27656)

# ------------------------------- Col-0 -------------------------------
print(col0.shape)                     # (2442, 27656)
print(issparse(col0.X))               # True
print(col0.X.dtype)                   # float64
print(list(col0.obs.columns))         # ['sample', 'n_genes_by_counts', 'total_counts',
                                      #  'total_counts_organellar', 'pct_counts_organellar', 'is_cell']
print(list(col0.var_names[:2]))       # ['AT1G01010.Araport11.447', 'AT1G01040.Araport11.447']
d = col0.X.data
d   = d[d != 0] 
print(f'{d.min():.4f} / {d.max():.4f}') ###0.0100 / 1075.8417

# Convert to int format, based on this ###0.0100 / 1075.8417, + the fact that we generated these counts 
# table these are raw counts
# convering them to int format.

# --- convert to integer counts ---------------------------------------
print('nnz before      :', col0.X.nnz)                    # 2287252

col0.X.data = np.round(col0.X.data)     # 0.01 -> 0,  2.7 -> 3
col0.X.eliminate_zeros()                # drop the values that rounded to 0
col0.X = col0.X.astype(np.int32)        # float64 -> int32

print('nnz after       :', col0.X.nnz)                    # fewer — the sub-0.5 dribs are gone
print('dtype           :', col0.X.dtype)                  # int32
print('min / max       :', col0.X.data.min(), col0.X.data.max())
print('all integers?   :', np.issubdtype(col0.X.dtype, np.integer))   # True

# --- save -------------------------------------------------------------
col0.write_h5ad(f'{MYB}/Col-0_INT_COUNTS.h5ad', compression='gzip')
print("^^^^^^^^^^^^Col-0 ARABIDOPSIS^^^^^^^^^^^^")

# ------------------------------- pFACT -------------------------------
print(pfact.shape)                    # (3304, 27656)
print(pfact.X.dtype)                  # float64

d_pfact = pfact.X.data
print(d_pfact.min())                  # 0.01
print(d_pfact.max())                  # 780.4822
print(np.allclose(d_pfact, np.round(d_pfact)))  # False
print(d_pfact.size)                   # 3682676
d = pfact.X.data[:]
d   = d[d != 0] 
print(f'{d.min():.4f} / {d.max():.4f}') ####0.0100 / 780.4822
# --- convert to integer counts ---------------------------------------
print('nnz before      :', pfact.X.nnz)                   # 3682676

pfact.X.data = np.round(pfact.X.data)    # 0.01 -> 0,  780.4822 -> 780
pfact.X.eliminate_zeros()                # prune the values that rounded to 0
pfact.X = pfact.X.astype(np.int32)       # float64 -> int32

print('nnz after       :', pfact.X.nnz)                   # fewer
print('dtype           :', pfact.X.dtype)                 # int32
print('min / max       :', pfact.X.data.min(), pfact.X.data.max())   # 1 780
print('integer dtype?  :', np.issubdtype(pfact.X.dtype, np.integer)) # True

# --- save -------------------------------------------------------------
pfact.write_h5ad(f'{MYB}/pFACT_INTCOUNTS.h5ad', compression='gzip')
print("^^^^^^^^^^^^pFACT ARABIDOPSIS NUCLEUS^^^^^^^^^^^^")

# ------------------------------ pHORST -------------------------------
print(phorst.shape)                   # (2937, 27656)
print(phorst.X.dtype)                 # float64

d_phorst = phorst.X.data
print(d_phorst.min())                 # 0.01
print(d_phorst.max())                 # 2415.6084
print(np.allclose(d_phorst, np.round(d_phorst)))  # False
print(d_phorst.size)                  # 4539900
d = phorst.X.data[:]
d   = d[d != 0] 
print(f'{d.min():.4f} / {d.max():.4f}') ###0.0100 / 2415.6084
# --- convert to integer counts ---------------------------------------
print('nnz before      :', phorst.X.nnz)                   # 4539900

phorst.X.data = np.round(phorst.X.data)   # 0.01 -> 0,  2415.6084 -> 2416
phorst.X.eliminate_zeros()                # prune the values that rounded to 0
phorst.X = phorst.X.astype(np.int32)      # float64 -> int32

print('nnz after       :', phorst.X.nnz)                   # fewer
print('dtype           :', phorst.X.dtype)                 # int32
print('min / max       :', phorst.X.data.min(), phorst.X.data.max())   # 1 2416
print('integer dtype?  :', np.issubdtype(phorst.X.dtype, np.integer))  # True

# --- save -------------------------------------------------------------
phorst.write_h5ad(f'{MYB}/pHORST_COUNTS.h5ad', compression='gzip')
print("^^^^^^^^^^^^pHORST ARABIDOPSIS NUCLEUS^^^^^^^^^^^^")

# ---------------------------- THE TRANSGENE --------------------------
# Sanity check that rounding did not eat the signal. These matrices are now int32,
# so the totals differ from the pre-rounding fractional values
# (Col-0 2.00 / pFACT 438.94 / pHORST 91.67). The CELL COUNTS are what must hold:
# 0 / 198 / 34. Col-0 going to 0 total is correct -- its 2.00 was entirely
# sub-1 fractional dribs, which is exactly why the detection threshold is >= 1 UMI.
i = col0.var_names.get_loc('AT4G28110.Fusion')
print(i)                              # column index, same in all three (shared reference)

tg_col0 = np.asarray(col0.X[:, i].todense()).ravel()
print(tg_col0.sum())                  # 0
print((tg_col0 >= 1).sum())           # 0     <- wild-type floor

tg_pfact = np.asarray(pfact.X[:, i].todense()).ravel()
print(tg_pfact.sum())
print((tg_pfact >= 1).sum())          # 198   <- must hold

tg_phorst = np.asarray(phorst.X[:, i].todense()).ravel()
print(tg_phorst.sum())
print((tg_phorst >= 1).sum())         # 34    <- must hold



