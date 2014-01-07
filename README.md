#Example

See mine/test_mine.py

```python
import numpy as np

x = np.array(range(1000))
y = 4 * (x - 1. / 2) ** 2
D = zip(x, y)
n = len(D)
B = pow(n, 0.6)
c = 15
M = ApproxCharacteristicMatrix(D, B, c=1)
print mine(M, B, c)
```
