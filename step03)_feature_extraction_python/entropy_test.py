import numpy as np
from scipy.stats import entropy
#p = [0.10,0.90,0.10,0.90,0.10,0.90,0.1,0.9,0.1,0.9,0.1,0.9,0.1,0.9]
#p = [0.5,0.15,0.5,0.85,0.9,0.5,0.9,0.95,0.95,0.85]
p=[0.8,0.9,0.8,0.8,0.7,0.7,0.2,0.2,0.2,0.2]
p = np.array(p)
normalized_entropy = entropy(p,base=2) / np.log2(p.shape[0])
print(normalized_entropy)