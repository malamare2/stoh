import numpy as np
import matplotlib.pyplot as plt

T8, C8, chi8 = np.loadtxt("f_T8.dat", unpack=True)
T16, C16, chi16 = np.loadtxt("f_T16.dat", unpack=True)
T32, C32, chi32 = np.loadtxt("f_T32.dat", unpack=True)

# graf
plt.scatter(T8, chi8, label="L = 8", s=15)
plt.scatter(T16, chi16, label="L = 16", s=15)
plt.scatter(T32, chi32, label="L = 32", s=15)

# osi
plt.xlabel("Temperatura T")
plt.ylabel("Magnetska susceptibilnost χ")

# NAĐI T gdje je maksimalna χ (uzima sve podatke zajedno)
T_all = np.concatenate([T8, T16, T32])
chi_all = np.concatenate([chi8, chi16, chi32])

idx_max = np.argmax(chi_all)
T_max = T_all[idx_max]

# naslov s T na maksimumu
plt.title(f"2D Ising model – susceptibilnost (T_max ≈ {T_max:.2f})")

# legenda i mreža
plt.legend()
plt.grid()



plt.show()