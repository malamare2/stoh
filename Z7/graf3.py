import numpy as np
import matplotlib.pyplot as plt

bm = []
Mm = []

with open("f_m8.dat", "r") as f:
    lines = f.readlines()

for line in lines[10003:20003]:

    if line.startswith("#"):
        continue

    p = line.split()

    bm.append(float(p[0]))
    Mm.append(float(p[1]))

bm = np.array(bm)
Mm = np.array(Mm)

Mcum = np.cumsum(Mm) / np.arange(1, len(Mm)+1)








# prorijedi podatke radi preglednosti
step = 5



# magnetizacija
plt.plot(
    bm[::step],
    Mm[::step],
    color='orange',
    lw=0.8,
    alpha=0.7,
    label=r'$M_b$ po bloku'
)

plt.plot(
    bm,
    Mcum,
    color='red',
    lw=2.5,
    label=r'$\langle M \rangle$ kumulativni prosjek'
)

plt.xlabel('Blok')
plt.ylabel('Energija')
plt.title('Magnetizacija tijekom simulacije')
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

