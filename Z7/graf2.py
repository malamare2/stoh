import numpy as np
import matplotlib.pyplot as plt

b = []
Eb = []

with open("f_E8.dat", "r") as f:
    lines = f.readlines()

# druga temperatura +uracunat burn-in
for line in lines[10003:20003]:

    if line.startswith("#"):
        continue

    p = line.split()

    #energija kao funkcija vremena simulacije
    b.append(float(p[0]))
    Eb.append(float(p[1]))

b = np.array(b)
#energija izmjenjena nakon svakog montecarlo bloka
Eb = np.array(Eb)

# kumulativni prosjek, zbaraja sve energije do tog trenutka 
#dijeli se brojem 
Ecum = np.cumsum(Eb) / np.arange(1, len(Eb)+1)

# prorijedi podatke radi preglednosti
step = 5

plt.figure(figsize=(12,5))

plt.plot(
    b[::step],
    Eb[::step],
    color='teal',
    lw=0.8,
    alpha=0.7,
    label=r'$E_b$ po bloku'
)

plt.plot(
    b,
    Ecum,
    color='darkblue',
    lw=2.5,
    label=r'$\langle E \rangle$ kumulativni prosjek'
)

plt.xlabel('Blok')
plt.ylabel('Energija')
plt.title('Energija tijekom simulacije')
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

