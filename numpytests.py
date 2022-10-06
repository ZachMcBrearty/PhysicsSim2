import numpy as np
np.zeros
a = np.zeros((5, 5))
b = np.zeros((5, 5))
f_a = np.vectorize(lambda x: x + 1)
f_b = np.vectorize(lambda x: 2 * (x + 3))

print(a)
print(f_a(a))
print(b)
print(f_b(b))