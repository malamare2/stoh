import numpy as np
import matplotlib.pyplot as plt

r, g = np.loadtxt("gr.dat", unpack=True)

plt.plot(r, g)
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radijalna distribucija")
plt.grid()
plt.show()