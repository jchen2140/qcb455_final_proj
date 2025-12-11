import numpy as np

data = np.load('GLS_p.npy')

# Print the first 10 elements
print(data[:10])

# Inspect useful metadata
print("\nShape of the array:", data.shape)
print("Data type (dtype):", data.dtype)
print("Number of elements:", data.size)