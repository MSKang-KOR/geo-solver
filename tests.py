
import numpy as np

a = [1,2,3,4,5]

def test(a):
    # b = np.array(a)
    # b = 10*b
    return np.array(a)

b = test(a)

print(type(a))
print(type(b))

a = [1,2]
print(a.index(-2) if -2 in a else a.index(2))
