import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(-4.1,4.1)
F=(np.exp(x)-1-x)/x**2

plt.plot(x,F)
plt.xlabel("$X$")
plt.ylabel("$ F(X) $")
plt.show()
