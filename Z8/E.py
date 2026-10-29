import numpy as np
import matplotlib.pyplot as plt

# Ucitavanje podataka o energiji
podaci = np.loadtxt("E.dat")

# Graf energije
plt.plot(podaci[:,0], podaci[:,1], label="Energija bloka")
plt.plot(podaci[:,0], podaci[:,2], label="Prosjek po blokovima")

plt.xlabel("Blok")
plt.ylabel("Energija")
plt.title("Energija tijekom simulacije")

plt.legend()
plt.grid()

plt.show()
