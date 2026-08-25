
import logging
import os
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

TRANSGENE = 'AT4G28110.Fusion'
ARAPORT_SUFFIX = '.Araport11.447'
RESOLUTION = 1.0
N_TOP_GENES = 3000
N_PCS = 30
SEED = 0

MYB = Path('/data1/bfernando/DE_project/alvinfry/h5ad')
COL0   = MYB / 'Col-0_INT_COUNTS.h5ad'
PFACT  = MYB / 'pFACT_INTCOUNTS.h5ad'
PHORST = MYB / 'pHORST_COUNTS.h5ad'

SAMPLES = {'Col-0': COL0, 
           'pFACT': PFACT, 
           'pHORST': PHORST}

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)-7s %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger('label')

def strip_ids(a: ad.AnnData) -> ad.AnnData:
    # isoforms are already collapsed to gene ID in alvinfry
    """AT1G01010.Araport11.447 -> AT1G01010, leaving the transgene untouched."""
    a.var_names = [g if g == TRANSGENE else g.replace(ARAPORT_SUFFIX, '')
                   for g in a.var_names]
    if a.var_names.has_duplicates:
        raise ValueError('ID strip produced duplicate gene names')
    if TRANSGENE not in a.var_names:
        raise KeyError(f'{TRANSGENE} lost during strip')
    return a

def cluster(a: ad.AnnData, name: str) -> ad.AnnData:
    """Leiden into dev_label, with .X returned to counts afterwards."""
    if not np.issubdtype(a.X.dtype, np.integer):
        raise ValueError(f'{name}: .X is {a.X.dtype}, expected integer counts')

    a.layers['_counts'] = a.X.copy()                  # park the counts

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_GENES)
    sc.tl.pca(a, n_comps=50, svd_solver='arpack', random_state=SEED)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=N_PCS, random_state=SEED)
    sc.tl.leiden(a, resolution=RESOLUTION, flavor='igraph', n_iterations=2,
                 random_state=SEED, key_added='dev_label')

    a.X = a.layers['_counts']                         # RESTORE -- do not skip
    del a.layers['_counts']
    a.uns.pop('log1p', None)                          # drop the stale stamp

    a.obs['dev_label'] = a.obs['dev_label'].astype(str).astype('category')
    log.info('%-7s %d clusters | .X back to %s | sizes %s', name,
             a.obs['dev_label'].nunique(), a.X.dtype,
             sorted(a.obs['dev_label'].value_counts().values, reverse=True))
    return a

def main() -> None:
    np.random.seed(SEED)
    parts = []

    for name, path in SAMPLES.items():
        a = ad.read_h5ad(path)
        a = strip_ids(a)
        a.obs['sample'] = name
        parts.append(a.copy())

        a = cluster(a, name)
        out = path.with_name(path.stem + '_LABELED.h5ad')
        a.write_h5ad(out, compression='gzip')
        log.info('wrote %s', out)

    # --- one-dataset layout: concatenate, then cluster once so the label space is shared
    joint = ad.concat(parts, label='sample', keys=list(SAMPLES), index_unique='-', join='inner')
    log.info('concatenated -> %s', joint.shape)
    joint = cluster(joint, 'joint')

    out = MYB / 'myb41_joint_LABELED.h5ad'
    joint.write_h5ad(out, compression='gzip')
    log.info('wrote %s', out)

    i = joint.var_names.get_loc(TRANSGENE)
    tg = np.asarray(joint.X[:, i].todense()).ravel()
    for s in SAMPLES:
        m = (joint.obs['sample'] == s).to_numpy()
        log.info('%-7s transgene >=1 UMI: %d', s, int((tg[m] >= 1).sum()))

