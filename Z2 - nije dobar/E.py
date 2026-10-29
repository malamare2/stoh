import numpy as np
import matplotlib.pyplot as plt
#slucajni hod jednak rezultat kao difuzijska jednadzba
# kako slučajno gibanje dovodi do difuzije

Nw = 64000
Nt = 200  
step = 3

L = 100   # velicina prostora od -100 do +100
dx = 2    # veličina kvadratića

nx = int(2*L/dx)   # koliko ima ćelija
ny = nx

x = np.zeros(Nw) #sve cestice krecu iz centra
y = np.zeros(Nw)

r2 = []  #prati se srednja kvadratna udaljenost tj koliko se sustav siri
steps = [50,100,199]  #trenutci u kojima spremamo raspodjelu
P_save = {}

for t in range(Nt):

    # slučajni pomak čestice svaka ga cestica dobiva u x, y smjeru a ravnomjerno gledamo jer svaki pomak ima jednaku vjerojatnost izmedu -3 i +3
    x += np.random.uniform(-step, step, Nw)
    y += np.random.uniform(-step, step, Nw)

    # srednja kvadratna udaljenost (udaljenost od centra)
    r2.append(np.mean(x**2 + y**2)) #to je prosjecna kvadratna udaljenost od centra
#teorija kaze r^2 = 4Dt

#e.1#test jeli random walk stvarno difuzija
#dobije se na posljetku linearni rast s vremenom


#sve cestice se mapiraju u diskretne celije mreze i nakon normalizacije imamo vjerojatnost
    # raspodjela P(x,y)
    #sve cest u mrezu celija i brojim koliko ih ima u celijama pretvaram vjerojatnost i crtam graf
    if t in steps:

        P = np.zeros((nx,ny))   #stvara mapu prostora

        for i in range(Nw):  #svaka cestica se gleda posebno

            # određivanje u kojoj je ćeliji čestica
            ix = int((x[i]+L)/dx)  #pretvara iz kontinuuma u diskretne mreze celija
            iy = int((y[i]+L)/dx) #pomak za L daje da nam budu poz cest a dx  odredujemo celiju kojoj pripada

            # provjera da je unutar prostora tj mreze
            if 0 <= ix < nx and 0 <= iy < ny:
                #distribucija
                 P[ix][iy] += 1   # broji koliko je čestica u ćeliji  ovo stvara diskretnu raspodjelu ili histogram

        P = P / Nw  #pretvara broj u vjerojatnost
        P_save[t] = P

t = np.arange(Nt)

plt.plot(t, r2)
plt.xlabel("t")
plt.ylabel("<r^2>")
plt.grid()
plt.show()


#matematicka provjera
#treba odrediti D
k,b = np.polyfit(t, r2, 1)  #nade pravac r2 =kt +b i provjerava linearnost
D = k/4  #ovo teorija kaze

print("D =", D)

# usporedba simulacije i teorijskog pravca
plt.plot(t, r2, label="simulacija")
plt.plot(t, k*t+b, "--", label="fit") #trazi linearnu aproksimaciju
plt.legend()
plt.grid()
plt.title("Fit <r^2> = 4Dt") #teorija kaze
plt.show()



plt.figure(figsize=(10,3))

#e.3 kako se raspodjela šetaća mijenja u vremenu

for i,t in enumerate(steps):
    plt.subplot(1,3,i+1) #dijeli sliku na mrezu 1 red 3 stupca
    plt.imshow(P_save[t], origin="lower")  #gdje su cestice u prostoru
    plt.title("t="+str(t))

plt.show()


#sad trazim difuzijsku jdz PDE
N = 60  #podjela prostora na mrezu kvadratica
rho = np.zeros((N,N))   #koliko ima u tom kvadraticu
rho[N//2][N//2] = 1   # početna gustoća u centru

D = 1
dt = 0.2

for t in range(200):  #200 trenutaka u vremenu simulira razvoj kroz vrijeme

    new = rho.copy() 

    for i in range(1,N-1): #ide kroz svaki kvadratic
        for j in range(1,N-1):  #da se izbjegnu rubovi

            # Laplace operator (utjecaj susjednih ćelija) gleda razliku centra i susjeda
            lap = (   #koliko je odstupa usporedba sa susjeduma i mjeri razliku okoline
                rho[i+1][j] + rho[i-1][j] + 
                rho[i][j+1] + rho[i][j-1] -
                4*rho[i][j]
            )

            # difuzijska jednadžba
            new[i][j] = rho[i][j] + D*dt*lap   #nova vrijednost u jednoj tocki mreze

    rho = new

#umerički rješava difuzijsku jednadžbu pomoću diskretizacije Laplaceovog operatora
plt.imshow(rho)
plt.title("PDE difuzija")
plt.colorbar()
plt.show()


plt.subplot(1,2,1)
plt.imshow(P_save[199])
plt.title("Random walk")

plt.subplot(1,2,2)
plt.imshow(rho)
plt.title("PDE")

plt.show()
#za kraj je usporedba random walka i PDE te daju slican rezultat
#usporedba je da pokazemo kajko mikroskopski model slucajni hod daje isti rezultat kao makroskopski model difuz jdz
