import numpy as np
import matplotlib.pyplot as plt

podaci = np.loadtxt("E_vs_a.dat")

a = podaci[:, 0]
E = podaci[:, 1]
sigma_E = podaci[:, 2]

indeks_min = np.argmin(E)
a_min = a[indeks_min]
E_min = E[indeks_min]
sigma_min = sigma_E[indeks_min]

Nw = 500

plt.figure(figsize=(11, 7))

plt.errorbar(
    a,
    E,
    yerr=100 * sigma_E,
    fmt="o-",
    markersize=5,
    linewidth=1.5,
    elinewidth=1.5,
    capsize=4,
    capthick=1.5,
    label=r"$E_{VMC}\;(\pm 100\sigma)$"
)

plt.title(
    fr"{Nw} šetača   "
    fr"(minimum: $a\approx{a_min:.2f}$, "
    fr"$E\approx{E_min:.2f}$)",
    fontsize=16
)

plt.xlabel(r"$a$", fontsize=14)
plt.ylabel(r"$E/\varepsilon_0$", fontsize=14)

plt.grid(alpha=0.3)
plt.legend(loc="upper left", fontsize=12)
plt.tick_params(axis="both", labelsize=12)

plt.tight_layout()
plt.savefig("energija_vs_a.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Optimalni parametar: a = {a_min:.2f}")
print(f"Najmanja energija: E = {E_min:.8f} ± {sigma_min:.8f}")