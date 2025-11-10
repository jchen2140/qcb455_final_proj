import numpy as np
import pandas as pd
import statsmodels.api as sm
from collections import defaultdict
from load_screens import load_screens
from scipy.stats import cauchy
from statsmodels.stats.multitest import fdrcorrection

screens = load_screens()
genes = screens.index

modules = []
d = 0.5
df = pd.read_csv(f"modules_d_{d}.csv", usecols=["Members"])

for members in df["Members"]:
    if pd.isna(members):
        continue
    gene_list = str(members).replace('"', '').split()
    if len(gene_list) > 0:
        modules.append(gene_list)

modules = np.array(modules, dtype=object)
print(f"Loaded {len(modules)} modules from modules_d_{d}.csv")

cell_line_cancer_types = (
    screens.columns.str.split("_").str[1:].str.join(" ").str.capitalize().to_series()
)
cancer_type_frequencies = cell_line_cancer_types.value_counts()
singletons = cancer_type_frequencies.index[cancer_type_frequencies == 1]
mask = ~cell_line_cancer_types.isin(singletons).values
screens = screens.iloc[:, mask]
cell_line_cancer_types = cell_line_cancer_types[mask]

# One-hot encode cancer types
one_hot_cancer_types = pd.get_dummies(cell_line_cancer_types).astype(float)
cancer_types = one_hot_cancer_types.columns
print(f"Analyzing {len(cancer_types)} cancer types")

cholsigmainv = np.linalg.cholesky(np.linalg.pinv(screens.cov())).T
inputs = cholsigmainv.dot(sm.add_constant(one_hot_cancer_types))

gene_ps = pd.DataFrame(index=genes, columns=cancer_types, dtype=float)

for gene in genes:
    output = cholsigmainv.dot(screens.loc[gene].astype(float))
    model = sm.OLS(output, inputs).fit()
    gene_ps.loc[gene] = model.pvalues[1:]

print("Computed gene-level p-values")

def ACAT(ps):
    ps = ps[~np.isnan(ps)]
    if len(ps) == 0:
        return np.nan
    return cauchy.sf(np.tan((0.5 - ps) * np.pi).mean())

results = defaultdict(list)

for gene_set in modules:
    gene_set_ps = gene_ps.loc[gene_ps.index.isin(gene_set)]
    for cancer_type in cancer_types:
        pvals = gene_set_ps[cancer_type].astype(float).values
        p_meta = ACAT(pvals)
        results["Module genes"].append(", ".join(gene_set))
        results["Cancer type"].append(cancer_type)
        results["p"].append(p_meta)

results = pd.DataFrame(results)
print(f"Computed p-values for {len(results)} module–cancer pairs")

results["FDR"] = results.groupby("Cancer type")["p"].transform(
    lambda p: fdrcorrection(p.fillna(1))[1]
)

filtered = results[results["FDR"] < 0.5].copy()
filtered = filtered.sort_values("p").reset_index(drop=True)
filtered.insert(0, "Rank", filtered.index + 1)

filtered.to_csv(
    "cancer_type_dependencies.tsv",
    sep="\t",
    index=False,
    float_format="%.2E",
    quotechar='"'
)
print(f"Saved {len(filtered)} moderately significant associations (FDR < 0.5)")
