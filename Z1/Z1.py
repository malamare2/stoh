import numpy as np
import matplotlib.pyplot as plt

#testira ima li korelacja imeđu brojeva udaljenoh za 3 mjesta

Nmax = 10000000
k = 3

x = np.random.random(Nmax + k)       #generira brojeve izmedu 0 i 1 i dobije se niz       (svi br izmedu su jednako vjerojatni)

koraci = []
test = []

zbroj = 0

for i in range(Nmax):                     #usporedujemo br koji su udaljeni za 3 mjesta

    zbroj += x[i] * x[i + k]   #dodavanje novog produkta
    C3 = zbroj / (i + 1)                    #trenutna vrijednost korelacije

    if (i + 1) % 10000 == 0:                #sluzi kako bi olaksala i prikazala manji br tocaka
        koraci.append(i + 1)
        test.append(np.sqrt(i + 1) * (C3 - 0.25))      #ako su sluv+cajni prosjek c3 treba bit 0.25

plt.plot(koraci, test)
plt.grid()
plt.xlabel("N")
plt.ylabel(r'$\sqrt{N}(C(3)-0.25)$')
plt.title("Susjed-susjed korelacije C(3)")
plt.show()


#prosjek produkta dva broja koja su nepovezama je 0.25 jer je prosjek da se dobije bilo koji izmedu 0 i 1 0.5
#Na malin N statistika je loša pa imaju velike fluktuacije a kako raste stabilizira se