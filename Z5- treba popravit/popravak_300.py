import numpy as np
import matplotlib.pyplot as plt

NbSkip = 50
Nb = 200
Nk = 100
Nw = 100
Nacc = 100

def Psi(r):
    return (27 - 18*r + 2*r**2) * np.exp(-r/3.)

ik = 0
ib = 1

x = np.zeros((4, Nw + 1))
xp = np.zeros(4)

dX = np.array([0, 2.0, 2.0, 2.0])

f = np.zeros(Nw + 1)
P = np.zeros(Nw + 1)

ff = open("r300.dat", "w")
frNw = open("rNw300.dat", "w")

# ---------------------------
# spremanje uzoraka (NOVO)
# ---------------------------
sample_x = []
sample_y = []
sample_r = []

# ---------------------------
# inicijalizacija
# ---------------------------
for iw in range(1, Nw + 1):
    rp2 = 0.0

    for k in range(1, 4):
        x[k][iw] = 20.0 * (np.random.rand() - 0.5)
        rp2 += x[k][iw] * x[k][iw]

    frNw.write("%7d  %13.8e  %13.8e  %13.8e\n" %
               ((ib - 1) * Nk + ik, x[1][iw], x[2][iw], x[3][iw]))

    rp = np.sqrt(rp2)
    P[iw] = Psi(rp)**2
    f[iw] = rp

frNw.write("\n\n")

# ---------------------------
# Monte Carlo
# ---------------------------
acc = 0.0
Sbf = 0.0
Sbf2 = 0.0

block_r = []
block_avg_list = []
sigma_blocks = []

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
                xp[k] = x[k][iw] + dx
                rp2 += xp[k] * xp[k]

            rp = np.sqrt(rp2)

            Pp = Psi(rp)**2
            T = Pp / P[iw]

            if T >= 1 or np.random.rand() <= T:
                for k in range(1, 4):
                    x[k][iw] = xp[k]
                acc += 1
                P[iw] = Pp
                f[iw] = rp

            Swf += f[iw]

        if ib > 0:
            Skf += Swf / Nw

    # ---------------------------
    # BLOCK ANALIZA
    # ---------------------------
    if ib > 0:

        block_avg = Skf / Nk

        Sbf += block_avg
        Sbf2 += block_avg**2

        ave_f = Sbf / ib
        sigma = np.sqrt(Sbf2 / ib - ave_f**2)

        sigma_blocks.append(sigma)
        block_r.append(block_avg)
        block_avg_list.append(ave_f)

        # ---------------------------
        # NOVO: spremanje uzoraka
        # ---------------------------
        if ib > Nb // 2:
            for iw in range(1, Nw + 1):
                sample_x.append(x[1][iw])
                sample_y.append(x[2][iw])
                sample_r.append(f[iw])

        ff.write("%7d  %13.8e  %13.8e\n" % (ib, block_avg, ave_f))

# ---------------------------
# konačni rezultati
# ---------------------------
ave_f = Sbf / Nb
sig_f = np.sqrt(Sbf2 / Nb - ave_f**2)
accP = acc / (Nb * Nw * Nk)

ff.close()
frNw.close()

print("\n accP = %4.1lf\n" % (accP * 100.0))
print("\n max dx = %6.4lf, %6.4lf, %6.4lf\n" % (dX[1], dX[2], dX[3]))
print("\n <r> = %8.5lf +- %8.5lf \n" % (ave_f, sig_f))

# ---------------------------
# GRAF 1: konvergencija
# ---------------------------
plt.figure()
plt.plot(block_r, label="block <r>")
plt.plot(block_avg_list, label="cumulative <r>")
plt.legend()
plt.xlabel("block")
plt.ylabel("<r>")
plt.title("Convergence of 3s orbital")
plt.grid()
plt.show()

# ---------------------------
# GRAF 2: fluktuacije
# ---------------------------
plt.figure()
plt.plot(block_r)
plt.xlabel("block")
plt.ylabel("<r>")
plt.title("Fluctuations")
plt.grid()
plt.show()

# ---------------------------
# GRAF 3: sigma blokovi
# ---------------------------
plt.figure()
plt.plot(sigma_blocks)
plt.xlabel("block")
plt.ylabel("sigma_r")
plt.title("Blocking error estimate")
plt.grid()
plt.show()

# ---------------------------
# NOVO GRAF 4: elektronski oblak
# ---------------------------
plt.figure(figsize=(6,6))
plt.scatter(sample_x, sample_y, s=1)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Electron cloud (3s orbital)")
plt.axis("equal")
plt.grid()
plt.show()

# ---------------------------
# NOVO GRAF 5: radijalna raspodjela
# ---------------------------
plt.figure()
plt.hist(sample_r, bins=50, density=True)
plt.xlabel("r")
plt.ylabel("Probability density")
plt.title(f"Radial distribution <r> = {ave_f:.3f} ± {sig_f:.3f}")
plt.grid()
plt.show()