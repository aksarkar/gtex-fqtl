import numpy as np
import pandas as pd
import tabix
import pyplink
import scipy.special as sp
import scipy.stats as st
import sys

_sf = st.norm().sf

def twas(twas_data, B, X):
  B = B.loc[twas_data['3_y']]
  # Hack: chr12 has a dupe
  B = B.loc[~B.index.duplicated(keep='first')]
  z = twas_data['5_x']
  z = z.loc[~z.index.duplicated(keep='first')]
  X_ = X.loc[:,twas_data['3_y']]
  X_ = X_.loc[~X_.index.duplicated(keep='first')]
  R = est_gwas_cov(X_)
  return z.dot(B) / np.sqrt(np.diag(B.T.dot(R).dot(B)))

def load_fqtl_model(prefix, pip_thresh=0.95):
  snp_annot = pd.read_csv(f'{prefix}/fqtl.snp.info.gz', sep='\t', header=None)
  tis_annot = pd.read_csv(f'{prefix}/fqtl.tis.info.gz', sep='\t', header=None)
  snp_lodds = pd.read_csv(f'{prefix}/fqtl.snp.lodds.txt.gz', sep='\t', header=None)
  keep = (sp.expit(snp_lodds) > pip_thresh).any(axis=0)
  L = pd.read_csv(f'{prefix}/fqtl.snp.theta.txt.gz', sep='\t', header=None)
  L.index = snp_annot[3]
  F = pd.read_csv(f'{prefix}/fqtl.tis.theta.txt.gz', sep='\t', header=None)
  F.index = tis_annot[1].apply(lambda x: x[:-9])
  B = L.loc[:,keep].dot(F.loc[:,keep].T)
  return snp_annot, B

def load_geno(prefix):
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
  return pd.DataFrame(vt.T.dot(np.diag(d)).dot(vt), index=x.columns, columns=x.columns)

if __name__ == '__main__':
  f = tabix.open(sys.argv[1])
  genes = pd.read_csv(sys.argv[2], sep='\t', header=None, index_col=0)
  result = dict()
  for k, row in genes.iterrows():
    try:
      snp_annot, B = load_fqtl_model(f'/broad/compbio/ypp/gtex/analysis/fqtl-gtex/result/fqtl-std/{k}/50/')
    except FileNotFoundError:
      print(f'warning: skipping empty model ({k}: {row[1]})')
      continue
    X = load_geno(f'/broad/hptmp/aksarkar/geno/{k}')
    try:
      gwas_z = (pd.DataFrame(f.query(f'chr{row[2]}', int(row[3] - 1e6), int(row[3] + 1e6)))
                .astype(dict(enumerate([str, int, int, str, str, float, float]))))
    except tabix.TabixError:
      print(f'warning: skipping empty cis-region in GWAS file {sys.argv[1]} ({k}: {row[1]})')
    try:
      twas_data = gwas_z.merge(snp_annot, left_on=[2, 3, 4], right_on=[3, 4, 5]).set_index('2_x')
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
    result[row[1]] = twas(twas_data, B, X)
  result = pd.DataFrame.from_dict(result, orient='index')
  result.to_csv(sys.argv[3], sep='\t', compression='gzip')
