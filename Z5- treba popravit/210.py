import numpy as np
import matplotlib.pyplot as plt

#parametri simulacije
NbSkip = 50  #pocetni blokovi koji se preskacu (nisu jos u ravnotezi)
Nb = 200  #br blokova koji se mjere
Nk = 100  #br koraka jednog bloka
Nw = 100   #br nasumičnih cestica tj walkera
Nacc = 100


def Psi(r, z):   #definira funkciju, r udaljenost od ishodista
    return z * np.exp(-r/2.)  #imamo eksponenc pad

#pocetnbe vrijednosti brojaca
ik = 0  #korak 
ib = 1  #blok


#matrica puna nula, 4 reda i 101 stupac
x = np.zeros((4, Nw + 1))
xp = np.zeros(4)  #privremena nova pozicija

#max pomak po koord
dX = np.array([0, 1.5, 1.5, 1.5])


f = np.zeros(Nw + 1)
P = np.zeros(Nw + 1)

#rezultati blokova
ff = open("r.dat", "w")
frNw = open("rNw.dat", "w")

#generiranje početnih položaja
for iw in range(1, Nw + 1):
    rp2 = 0.0
    for k in range(1, 4):

        #svaki walker ima slučajnu poč poziciju pa imamo mnogo različčitih poč položaja
        x[k][iw] = 14.0 * (np.random.rand() - 0.5)
        rp2 += x[k][iw] * x[k][iw]

    frNw.write("%7d  %13.8e  %13.8e  %13.8e\n" %
               ((ib - 1) * Nk + ik, x[1][iw], x[2][iw], x[3][iw]))

    #racun vjerojatnosti za svakog setaca, stara vjerojstnost
    rp = np.sqrt(rp2)
    tmp = Psi(rp, x[3][iw])
    P[iw] = tmp * tmp
    f[iw] = rp

frNw.write("\n\n")

#monte carlo
acc = 0.0
Sbf = 0.0
Sbf2 = 0.0

block_r = []
block_avg_list = []
sigma_blocks = []  

#krece od negat bloka da se stigne zagrijat tj da mozr prve odbaciti
for ib in range(1 - NbSkip, Nb + 1):

    Skf = 0.0

    if ib == 1:
        acc = 0.0

    for ik in range(1, Nk + 1):

        Swf = 0.0

        for iw in range(1, Nw + 1):

            rp2 = 0.0

            for k in range(1, 4):
                dx = 2.0 * (np.random.rand() - 0.5) * dX[k]
                
                #metropolis korak, nova pozicija
                xp[k] = x[k][iw] + dx
                rp2 += xp[k] * xp[k]
            
            #nova vjerojatnost
            rp = np.sqrt(rp2)
            tmp = Psi(rp, xp[3])
            #racunanje nove vjerojatnosti
            Pp = tmp * tmp

            #metropolis kriterij
            T = Pp / P[iw]

            #ako je nova pozicija vjerojatnija prihvati se uvjet
            #ako je lošija ponekad prihvati
            if T >= 1 or np.random.rand() <= T:
                for k in range(1, 4):
                    x[k][iw] = xp[k]
                acc += 1
                P[iw] = Pp
                #mjerenje udaljenosti elektrona od jezgre
                #svaki walker ima udaljenost r
                f[iw] = rp

            #suma svih walkera
            Swf += f[iw]

        #prosjek po bloku
        if ib > 0: #brojac blokova
            Skf += Swf / Nw  #suma svih udaljenosti r, prosjek kroz cijeli blok


    if ib < 0:  #prvi blokovi simulacije  (1-NbSkip)

        #koliko cesto se prihvacaju novi polozaji       
        accP = acc / (Nk * Nw)
        acc = 0.0

        #povecava i smanjuje korake obzirom koliko prihvaca
        if accP > 0.6:
            dX *= 1.05
        elif accP < 0.4:
            dX *= 0.95


# BLOCK ANALIZA
    if ib > 0:
#prosjek bloka
        block_avg = Skf / Nk

        #srednja vrijednost i kvadrat srednje vrijednosti
        Sbf += block_avg
        Sbf2 += block_avg**2

        ave_f = Sbf / ib

        #statistička nesigurnost rezultata(fluktacije)
        sigma = np.sqrt(Sbf2 / ib - ave_f**2)
        sigma_blocks.append(sigma)

        block_r.append(block_avg)
        block_avg_list.append(ave_f)

        if ib % (Nb // 10) == 0:
            print("ib = %d" % ib)

        ff.write("%7d  %13.8e  %13.8e\n" % (ib, block_avg, ave_f))

#prosjek svih blokova do sad
ave_f = Sbf / Nb
sig_f = np.sqrt(Sbf2 / Nb - ave_f**2)
accP = acc / (Nb * Nw * Nk)

ff.close()
frNw.close()

print("\n accP = %4.1lf\n" % (accP * 100.0))
print("\n max dx = %6.4lf, %6.4lf, %6.4lf\n" % (dX[1], dX[2], dX[3]))
print("\n <r> = %8.5lf +- %8.5lf \n" % (ave_f, sig_f))

plt.figure()
plt.plot(block_r, label="block <r>")
plt.plot(block_avg_list, label="cumulative <r>")
plt.axhline(5.0, linestyle="--", label="exact (5)")
plt.legend()
plt.xlabel("block")
plt.ylabel("<r>")
plt.title("Convergence of <r>")
plt.grid()
plt.show()

plt.figure()
plt.plot(block_r)
plt.xlabel("block")
plt.ylabel("<r> per block")
plt.title("Block fluctuations")
plt.grid()
plt.show()

plt.figure()
plt.plot(sigma_blocks)
plt.xlabel("block")
plt.ylabel("sigma_r")
plt.title("Standard deviation (blocking method)")
plt.grid()
plt.show()


