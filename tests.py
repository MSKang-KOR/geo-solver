
import numpy as np

a = np.ones(5).reshape(5,1)
b = np.zeros(5).reshape(5,1)

# print(np.shape(a))
# print(np.shape(a.T))
c = np.concatenate((a,b),axis=1)
print(c)
