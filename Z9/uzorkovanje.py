import numpy as np
import matplotlib.pyplot as plt

podaci = np.loadtxt("Energija.dat")

blok = podaci[:, 0]
E_blok = podaci[:, 1]
E = podaci[:, 2]

E0 = E[-1]
sigma = np.std(E_blok, ddof=1) / np.sqrt(len(E_blok))

alpha = 0.6      # optimalni alpha
Nw = 500
Nt = 200

plt.figure(figsize=(11,6))

plt.plot(blok, E_blok,
         color="cyan",
         lw=1,
         label=r"$\langle E_b\rangle$")

plt.plot(blok, E,
         color="blue",
         lw=2.5,
         label=r"$\langle E\rangle$")

plt.title(
    fr"{Nw} šetača, VMC ($\alpha={alpha:.2f}$)    "
    fr"$E_0={E0:.3f}\pm{sigma:.3f}$"
)

plt.xlabel(f"blokovi / {Nt} koraka")
plt.ylabel(r"$E/\varepsilon_0$")

plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("energija_ravnotezno_uzorkovanje.png", dpi=300)
plt.show()