# Written by Yubi

import numpy as np
from load_screens import load_screens
from scipy.special import stdtr

# Load batch-corrected screens

screens = load_screens()

# Remove cell lines with any missing genes
# (not required for DepMap 18Q3, but is for more recent releases)
# You can use other strategies to remove NaNs instead, like imputing,
# removing genes with any missing cell lines

screens.dropna(axis=1, inplace=True)

# save screens as a csv file
screens.to_csv("bias_corrected_data.csv")
