r"""TWAS using fqtl models

The major differences between the implementation here and existing
implementations in S-PREDIXCAN (Barbeira et al. 2018) and FUSION (Gusev et
al. 2016, Mancuso et al. 2017) are:

1. We estimate the TWAS statistic based on the following argument for tissue t:

   y \sim N(eta a, sigma^2 I),    eta = X b
   \hat{a} = (eta' eta)^{-1} eta' y
           = (b' X' X b)^{-1} b' X' y
           = (b' R b)^{-1} b' z,    R = X' X / n, z = X' y / n

2. We hard-threshold tissue effect sizes based on PIP, in order to avoid
   correlating small posterior mean effect sizes with positive GWAS signals.

2. We use adaptive permutations (Che et al. 2014) to estimate p-values, rather
   than making a Normal assumption under the null.

"""
import numpy as np
import pandas as pd
import tabix
import pyplink
import scipy.sparse as ss
import scipy.special as sp
import scipy.stats as st
import sys

_sf = st.norm().sf

def align(twas_data, X, B, pip):
  """Align the TWAS data

  Return:

  z - GWAS z-scores (p,)  [Important: this needs to be broadcastable]
  B - fqtl effect size matrix (p, m)
  pip - fqtl posterior inclusion probabilities (p, m)

  """
  # Align the TWAS data, genotypes, and fqtl effect size matrix
  B = B.loc[twas_data['3_y']]
  # Hack: chr12 has a dupe in fqtl fitted models "12_48000000_T_C_b37: 2 times"
  B = B.loc[~B.index.duplicated(keep='first')]
  pip = pip.loc[twas_data['3_y']]
  pip = pip.loc[~pip.index.duplicated(keep='first')]
  z = twas_data['5_x']
  z = z.loc[~z.index.duplicated(keep='first')]
  X_ = X.loc[:,~X.columns.duplicated(keep='first')]
  X_ = X_.loc[:,twas_data['3_y']]
  X_ = X_.loc[:,~X_.columns.duplicated(keep='first')]
  R = est_gwas_cov(X_)
  return z, R, B, pip

def adaptive_permutation(z, R, B, pip, pip_thresh=0.95, max_success=10, batch_size=1000, max_perm=10_000_000):
  """Return TWAS approximate permutation p-values

  Return:

  pvals - p-values (m, 1)

  """
  B = B[pip > pip_thresh].fillna(0)
  n_qtls = (pip > pip_thresh).sum(axis=0)
  num_success = 0
  num_perm = 0
  stats = []
  pvals = []
  for t in B:
    if n_qtls[t] == 0:
      stats.append(float(0))
      pvals.append(float(1))
      continue
    stat = (z @ B[t]) / (B[t].T @ R @ B[t])
    stats.append(stat)
    # Generate permuted TWAS statistics in batches by using sparse matrices to
    # represent a batch of permuted non-zero effects
    for _ in range(max_perm // batch_size):
      perm_stat = []
      data = np.repeat(B[t][pip[t] > pip_thresh], batch_size)
      row = np.random.choice(B.shape[0], replace=True, size=(n_qtls[t], batch_size)).ravel()
      col = np.repeat(np.arange(batch_size), n_qtls[t])
      bt = ss.coo_matrix((data, (row, col)), shape=(B.shape[0], batch_size)).tocsr()
      num = (z.values.reshape(1, -1) @ bt)
      den = (R.values @ bt.power(2)).sum(axis=0, keepdims=True)
      perm_stat.append(num / den)
      perm_stat = np.vstack(perm_stat)

      num_success += (np.square(perm_stat) >= np.square(stat)).sum()
      num_perm += batch_size
      if num_success >= max_success:
        pvals.append(num_success / num_perm)
        break
    else:
      pvals.append((num_success + 1) / (max_perm + 1))
  # Important: after filtering low PIP factors/loadings, B is not guaranteed to
  # have all 44 tissues
  stats = pd.Series(stats, index=B.columns)
  pvals = pd.Series(pvals, index=B.columns)
  return stats, pvals

def load_fqtl_model(prefix, pip_thresh=0.95):
  """Return fqtl metadata, effect sizes, and PIPs

  Return:

  snp_annot - DataFrame (chromosome, ID, dummy, position, ref, alt)
  B - fqtl effect sizes (p, m)
  pip - fqtl posterior inclusion probabilities (p, m)

  """
  snp_annot = pd.read_csv(f'{prefix}/fqtl.snp.info.gz', sep='\t', header=None)
  snp_lodds = pd.read_csv(f'{prefix}/fqtl.snp.lodds.txt.gz', sep='\t', header=None)
  snp_lodds.index = snp_annot[3]
  tis_annot = pd.read_csv(f'{prefix}/fqtl.tis.info.gz', sep='\t', header=None)
  tis_lodds = pd.read_csv(f'{prefix}/fqtl.tis.lodds.txt.gz', sep='\t', header=None)
  tis_lodds.index = tis_annot[1].apply(lambda x: x.rsplit('_', maxsplit=1)[0])
  keep = (sp.expit(snp_lodds) > pip_thresh).any(axis=0)
  L = pd.read_csv(f'{prefix}/fqtl.snp.theta.txt.gz', sep='\t', header=None)
  L.index = snp_annot[3]
  F = pd.read_csv(f'{prefix}/fqtl.tis.theta.txt.gz', sep='\t', header=None)
  F.index = tis_annot[1].apply(lambda x: x[:-9])
  B = L.loc[:,keep].dot(F.loc[:,keep].T)
  pip = pd.DataFrame(sp.expit(snp_lodds.loc[:,keep]) @ sp.expit(tis_lodds.loc[:,keep]).T,
                     index=B.index, columns=B.columns)
  assert np.isfinite(pip.values).all()
  return snp_annot, B, pip

def load_sqtl_model(prefix, pip_thresh=0.95):
  """Return sqtl metadata, effect sizes, and PIPs

  Return:

  snp_annot - DataFrame (chromosome, ID, dummy, position, ref, alt)
  B - fqtl effect sizes (p, m)
  pip - fqtl posterior inclusion probabilities (p, m)

  """
  tis_annot = pd.read_csv(f'{prefix}/sqtl.tis.info.gz', sep='\t', header=None, index_col=0)
  snp_annot = pd.read_csv(f'{prefix}/sqtl.snp.info.gz', sep='\t', header=None)
  snp_lodds = pd.read_csv(f'{prefix}/sqtl.effect.lodds.txt.gz', sep='\t', header=None)
  pip = sp.expit(snp_lodds)
  pip.index = snp_annot[3]
  pip.columns = tis_annot[1].apply(lambda x: x.rsplit('_', maxsplit=1)[0])
  B = pd.read_csv(f'{prefix}/sqtl.effect.theta.txt.gz', sep='\t', header=None)
  B.index = snp_annot[3]
  B.columns = tis_annot[1].apply(lambda x: x.rsplit('_', maxsplit=1)[0])
  return snp_annot, B, pip

def load_geno(prefix):
  """Return DataFrame containing genotypes"""
  with pyplink.PyPlink(prefix) as f:
    fam = f.get_fam()
    bim = f.get_bim()
    x = np.ma.masked_equal(np.array([row for _, row in f.iter_geno()]).T, -1).astype(float)
    x -= x.mean()
    x /= x.std()
    x = pd.DataFrame(x, index=fam['fid'], columns=bim['pos'])
    return x

def est_gwas_cov(x, sv_thresh=1e-4):
  u, d, vt = np.linalg.svd(x, full_matrices=False)
  d[d < sv_thresh] = 0
  return pd.DataFrame(vt.T.dot(np.diag(np.square(d))).dot(vt) / x.shape[0], index=x.columns, columns=x.columns)

if __name__ == '__main__':
  f = tabix.open(sys.argv[1])
  genes = pd.read_csv(sys.argv[2], sep='\t', header=None, index_col=0)
  fqtl_stats = dict()
  fqtl_pvals = dict()
  sqtl_stats = dict()
  sqtl_pvals = dict()
  for k, row in genes.iterrows():
    try:
      snp_annot, B, pip = load_fqtl_model(f'/broad/compbio/ypp/gtex/analysis/fqtl-gtex/result/fqtl-std/{k}/50/')
    except FileNotFoundError:
      print(f'warning: skipping empty model ({k}: {row[1]})')
      continue

    try:
      gwas_z = (pd.DataFrame(f.query(f'chr{row[2]}', int(row[3] - 1e6), int(row[3] + 1e6)))
                .astype(dict(enumerate([str, int, int, str, str, float, float]))))
    except tabix.TabixError:
      print(f'warning: skipping empty cis-region in GWAS file {sys.argv[1]} ({k}: {row[1]})')
      continue

    try:
      twas_data = gwas_z.merge(snp_annot, left_on=2, right_on=3).set_index(2)
      strand_flip = dict(zip('ACGT', 'TGCA'))
      temp = twas_data[['3_x', '4_x']].applymap(lambda x: strand_flip[x])
      strand_flipped = np.logical_or(
        np.logical_and(temp['3_x'] == twas_data['4_y'], temp['4_x'] == twas_data['5_y']),
        np.logical_and(temp['3_x'] == twas_data['5_y'], temp['4_x'] == twas_data['4_y']))
      twas_data.loc[strand_flipped, ['3_x', '4_x']] = temp
      allele_flipped = np.logical_and(twas_data['3_x'] == twas_data['5_y'], twas_data['4_x'] == twas_data['4_y'])
      twas_data.loc[allele_flipped, '3_x'], twas_data.loc[allele_flipped, '4_x'] = \
        twas_data.loc[allele_flipped, '4_x'], twas_data.loc[allele_flipped, '3_x']
      twas_data.loc[allele_flipped, '5_x'] *= -1
      aligned = np.logical_and(twas_data['3_x'] == twas_data['4_y'], twas_data['4_x'] == twas_data['5_y'])
      twas_data = twas_data.loc[aligned]
    except:
      # Hack: genes like ENSG00000215784.4 have garbage in snp_annot
      print(f'warning: skipping corrupted model ({k}: {row[1]})')
      continue

    if twas_data.empty:
      print(f'warning: no SNPs left ({k}: {row[1]})')
      continue
    elif twas_data.shape[0] < 0.5 * snp_annot.shape[0]:
      print(f'warning: many SNPs lost ({k}: {row[1]})')
    else:
      print(f'estimating sTWAS statistics ({k}: {row[1]})')
      X = load_geno(f'/broad/hptmp/aksarkar/geno/{k}')
      z, R, B, pip = align(twas_data, X, B, pip)
      stat, pval = adaptive_permutation(z, R, B, pip)
      fqtl_stats[row[1]] = stat
      fqtl_pvals[row[1]] = pval

    try:
      snp_annot, B, pip = load_sqtl_model(f'/broad/compbio/ypp/gtex/analysis/fqtl-gtex/result/sqtl/{row[2]}/{row.name}/')
    except FileNotFoundError:
      print(f'warning: skipping empty sqtl model ({k}: {row[1]})')
      continue
    z, R, B, pip = align(twas_data, X, B, pip)
    stat, pval = adaptive_permutation(z, R, B, pip)
    sqtl_stats[row[1]] = stat
    sqtl_pvals[row[1]] = pval

  pd.DataFrame.from_dict(sqtl_stats, orient='index').to_csv(f'{sys.argv[3]}.sqtl.stat.txt.gz', sep='\t', compression='gzip')
  pd.DataFrame.from_dict(sqtl_pvals, orient='index').to_csv(f'{sys.argv[3]}.sqtl.pval.txt.gz', sep='\t', compression='gzip')
  pd.DataFrame.from_dict(fqtl_stats, orient='index').to_csv(f'{sys.argv[3]}.fqtl.stat.txt.gz', sep='\t', compression='gzip')
  pd.DataFrame.from_dict(fqtl_pvals, orient='index').to_csv(f'{sys.argv[3]}.fqtl.pval.txt.gz', sep='\t', compression='gzip')
