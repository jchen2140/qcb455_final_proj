import numpy as np
import pandas as pd
import statsmodels.api as sm
from collections import defaultdict
from load_screens import load_screens
from scipy.stats import cauchy
from statsmodels.stats.multitest import fdrcorrection

# Load screens
screens = load_screens()
genes = screens.index

# Load modules properly
modules = []
for d in (0.2, 0.5, 0.9):
    df = pd.read_csv(f'modules_d_{d}.csv', usecols=['Members'])
    for members in df['Members']:
        gene_list = str(members).replace('"', '').split()
        modules.append(gene_list)
modules = np.array(modules, dtype=object)

# Get cancer types; filter out singletons
cell_line_cancer_types = screens.columns.str.split('_').str[1:].str.join(' ').str.capitalize().to_series()
cancer_type_frequencies = cell_line_cancer_types.value_counts()
singletons = cancer_type_frequencies.index[cancer_type_frequencies == 1]
mask = ~cell_line_cancer_types.isin(singletons).values
screens = screens.iloc[:, mask]
cell_line_cancer_types = cell_line_cancer_types[mask]

# One-hot encode
one_hot_cancer_types = pd.get_dummies(cell_line_cancer_types).astype(float)
cancer_types = one_hot_cancer_types.columns

# GLS inputs
cholsigmainv = np.linalg.cholesky(np.linalg.pinv(screens.cov())).T
inputs = cholsigmainv.dot(sm.add_constant(one_hot_cancer_types))

# Compute gene-level p-values
gene_ps = pd.DataFrame(index=genes, columns=cancer_types, dtype=float)
for gene in genes:
    output = cholsigmainv.dot(screens.loc[gene].astype(float))
    model = sm.OLS(output, inputs).fit()
    gene_ps.loc[gene] = model.pvalues[1:]

# ACAT function
ACAT = lambda ps: cauchy.sf(np.tan((0.5 - ps) * np.pi).mean())

# Compute module-level p-values
results = defaultdict(list)
for gene_set in modules:
    gene_set_ps = gene_ps.loc[gene_ps.index.isin(gene_set)]
    for cancer_type in cancer_types:
        p_meta = ACAT(gene_set_ps[cancer_type].astype(np.float64))
        results['Module genes'].append(', '.join(gene_set))
        results['Cancer type'].append(cancer_type)
        results['p'].append(p_meta)

results = pd.DataFrame(results)

# FDR correction by cancer type
results['FDR'] = results.groupby('Cancer type')['p'].transform(lambda p: fdrcorrection(p)[1])

# Rank and save
results = results.sort_values('p').reset_index(drop=True)
results.insert(0, 'Rank', results.index + 1)

results.to_csv('cancer_type_dependencies.tsv', sep='\t', index=False, float_format='%.2E', quotechar='"')
