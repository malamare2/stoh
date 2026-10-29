import numpy as np
import matplotlib.pyplot as plt


dx = dy = 2 #podjela prostora na kvadratice

Nw = 800000 
Nk = 4000 #broj koraka (daje razvoj u vremenu)

D = 1e-2 #difuzijski koef
dt = 1 #jedan vremenski korak
kB = 1 #boltzmannova konst

x = np.random.uniform(40, 60, Nw) #po;etne poyicije malog kvadrata di se nalaze cestice
y = np.random.uniform(40, 60, Nw)

nx = int(100 / dx)  #broj celija, prostor je 100*100
ny = int(100 / dy)

save_steps = [50, 200, 500, 1500]  #trenutci koje spremamo tj broj vremenskih koraka k
results = {}

entropy = []
time = []

#petlja koja vodi evoluciju sustava
for k in range(Nk):     #svaki prolay je jednak jednom vremenskom koraku


    # Brownovo gibanje - nije ravnomjerno kretanje nego puno sitnih udara koji daju normalnu gausovu raspodjelu pomaka
    #0 imamo zapisan jer nema preferirani smjer cestica ide lijevo desno gore i dolje jednako vjerojkatno
    x += np.random.normal(0, np.sqrt(2*D*dt), Nw) #svaka cest dobije sslucajni pomak neovisno o smjeru
    y += np.random.normal(0, np.sqrt(2*D*dt), Nw) #devijacija
    #srednje kvadratno pomicanje raste linearno s vremenom

    # refleksija na rubovima
    x[x < 0] = -x[x < 0]
    x[x > 100] = 200 - x[x > 100]

    y[y < 0] = -y[y < 0]
    y[y > 100] = 200 - y[y > 100]

#prostor se dijeli na celije i brojim cestice u svakoj celiji
    # histogram (gustoća)   distribucija
    P, _, _ = np.histogram2d( #Pij je broj cestica u celiji
        x, y, 
        bins=[nx, ny], #nx x ny je prostor podjeljen na mrežu
        range=[[0, 100], [0, 100]]
    )

#ovaj racuna vjerojatnost i iz toga entropiju  tj neuredenje 
    P = P / Nw #vjerojatnost da je cestica u celiji   vjerojatnost potrebnu za entropiju

    # entropija
    Pnz = P[P > 0]  #uzima celije di ima cestica
    S = -kB * np.sum(Pnz * np.log(Pnz)) #shannon boltzmanova entopija
    #Pnz je vjerojatnost cestice u celiji a lnP mjera rijetkosti

    entropy.append(S)
    time.append(k)

    # spremanje za prikaz
    if k in save_steps:
        results[k] = P


#grid koordinata za graf
x_cent = np.linspace(0, 100, nx)
y_cent = np.linspace(0, 100, ny)
X, Y = np.meshgrid(x_cent, y_cent)



fig = plt.figure(figsize=(10, 7))

#crtam 4 grafa razlicitih vremena
for i, k in enumerate(save_steps):

    ax = fig.add_subplot(2, 2, i+1, projection='3d')

    Z = results[k].T  

    ax.plot_surface(X, Y, Z, cmap="Reds")

    ax.set_title("k = " + str(k))
    ax.set_xlabel("x")
    ax.set_ylabel("y")

plt.tight_layout()
plt.show()


#entropija
plt.figure()
plt.plot(time, entropy)
plt.xlabel("k")
plt.ylabel("S")
plt.title("Entropija sustava")
plt.grid()
plt.show()