import numpy as np
import matplotlib.pyplot as plt


# mikrokanonski sustav
N = 100                      # broj čestica
ES = 20                      # ukupna energija sustava
nstep = 1000000              # broj Monte Carlo koraka
nskip = 100000               # broj koraka za stabilizaciju
dvMax = 1.0                  # maksimalna promjena brzine


# početne brzine
v0 = np.sqrt(2 * ES / N)
#lista  brzina, svaka cet ima istu poc brzinu v0
v = []
for i in range(N):
    v.append(v0)
v = np.array(v)

ED = 0.0                     # energija demona
ED_list = []                 # lista energija demona kroz vrijeme


# kinetička energija
def E(v):
    return 0.5 * v**2 #što je veća v veća je i E






# Monte Carlo demon algoritam
for step in range(nstep):
    # birsmo jednu nasumično odabrana čestica
    ip = np.random.randint(N)
    # slučajna promjena brzine iz [-dvMax, dvMax]
    # dakle v se moze povecati ili smanjiti
    dv = dvMax * (2 * np.random.rand() - 1)
    # nova brzina
    vTmp = v[ip] + dv
    # promjena energije
    dE = E(vTmp) - E(v[ip])

    # demon prihvaća ili odbija promjenu
    #ako ima dovoljno energ moze je platiti
    if ED >= dE:
        ED = ED - dE
        #ako je uvjet ispunjen prihvaca brzinu
        v[ip] = vTmp

    # nakon stabilizacije spremamo energiju demona
    if step >= nskip:
        ED_list.append(ED)


ED_list = np.array(ED_list)
# ukupna energija je očuvana
ES_list = ES - ED_list



#PRIPREMA PODATAKA

# x-os, broji korake 1,2,3... tj stvara listu zacrtanje grafa
##promjena sustava kroz vrijeme
istep = [] #stvaram listu jer nam treba x-os (vrijeme)
for i in range(len(ED_list)):
    istep.append(i + 1)

# tekući prosjeci za  demona
ED_avg_run = []
suma = 0

#racunbamoi prosjek sa smanjnenim šumom
#svakom koraku se dodaje nova vrijednost demona
for i in range(len(ED_list)):
    suma += ED_list[i]   #dodavanje nove vrijednosti energije demona
    ED_avg_run.append(suma/(i+1)) #prosjek svih do sad uracunatih
    

#za sustav prosjek energije
ES_avg_run = []
suma = 0

for i in range(len(ES_list)):
    suma += ES_list[i]
    ES_avg_run.append(suma/(i+1))






#energija sustava
plt.figure()
plt.plot(istep, ES_list, alpha=0.4)
plt.plot(istep, ES_avg_run, 'r')
plt.title("Energija sustava")
plt.xlabel("koraci")
plt.ylabel("ES")
plt.show()


#energija demona
plt.figure()
plt.plot(istep, ED_list, alpha=0.4)
plt.plot(istep, ED_avg_run, 'r')
plt.title("Energija demona")
plt.xlabel("koraci")
plt.ylabel("ED")
plt.show()
#sustav skace ali je prosjek stabilamn


#zbrajanje svih vrijednosti, prosjek
suma = 0
for i in range(len(ED_list)):
    suma += ED_list[i]

#prosjhek energije demona
ED_avg = suma / len(ED_list)
ES_avg = ES - ED_avg

#spremanje rezulatat
tablica = {
    "N": N,
    "ES": ES,
    "<ED>": ED_avg,
    "<ES>": ES_avg,
    "<ES>/N": ES_avg / N,
    "0.5*N*<ED>": 0.5 * N * ED_avg #teorijaks provjera
}



with open("results.txt", "w") as f:
    for key, value in tablica.items():
        f.write(f"{key:15s}{value}\n")

#raspodjela energije, koliko cesto ima demon odredenu energiju, vjerojatnost
hist, bins = np.histogram(ED_list, bins=100, density=True)
ED_center = (bins[:-1] + bins[1:]) / 2


plt.figure()
plt.plot(ED_center, hist, 'r')
plt.title("P(ED)")
plt.xlabel("ED")
plt.ylabel("P(ED)")
plt.show()

#imam sve vrijednosti demona i djelim ih u intervale
#koliko cesto se ED pojavljuje u svakom intervalu
hist, bins = np.histogram(ED_list, bins=30, density=True)
#daje tocke u sredini svakog intervala za crtat graf
ED_center = (bins[:-1] + bins[1:]) / 2

#izbacuje lose vrijednosti, ne moze biti log(0)
mask = (hist > 0) & (ED_center <= 3)

#teorija kaze P = exp(-ED/T), pa log daje pravac y = mx +b
x = ED_center[mask]  #energija demona
y = np.log(hist[mask]) # log vjerojatnost


# linearni fit, pronalazi se nagib, najbolji pravac
m_fit, b_fit = np.polyfit(x, y, 1)


plt.figure()
#visoki T sporo padanje i obrnuto
#sto je veca energija manja je vjerojatnost pa je pad ekspon.
plt.scatter(x, y, s=10)
plt.plot(x, m_fit * x + b_fit)
plt.title("ln P(ED)")
plt.xlabel("ED")
plt.ylabel("ln P(ED)")
plt.grid()
plt.show()


#temperatura
print("nagib =", m_fit) 
print("T =", -1 / m_fit)