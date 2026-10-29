import numpy as np    
import matplotlib.pyplot as plt


dx = dy = 2  #podjela prostora na kvadratice

Nw = 800000 
Nk = 4000 #broj koraka (daje razvoj u vremenu)

D = 1e-2   #difuzijski koef, određuje jačinu slučajnih pomaka čestica
dt = 1   #jedan vremenski korak


x = np.random.uniform(40, 60, Nw)   #po;etne poyicije malog kvadrata di se nalaze cestice
y = np.random.uniform(40, 60, Nw)

nx = int(100 / dx)     #broj celija, prostor je 100*100
ny = int(100 / dy)

save_steps = [50, 200, 500, 1500]  #trenutci koje spremamo tj broj vremenskih koraka k
results = {}


#petlja koja vodi evoluciju sustava
for k in range(Nk):     #svaki prolay je jednak jednom vremenskom koraku

    # Brownovo gibanje - nije ravnomjerno kretanje nego puno sitnih udara koji daju normalnu gausovu raspodjelu pomaka
    #0 imamo zapisan jer nema preferirani smjer cestica ide lijevo desno gore i dolje jednako vjerojkatno
    x += np.random.normal(0, np.sqrt(2*D*dt), Nw)     #svaka cest dobije sslucajni pomak neovisno o smjeru
    y += np.random.normal(0, np.sqrt(2*D*dt), Nw) 
    #srednje kvadratno pomicanje raste linearno s vremenom 


    # refleksija na rubovima ako cest izade iz kutije odbije se neazad pa nema gubitka cestica
    x[x < 0] = -x[x < 0]
    x[x > 100] = 200 - x[x > 100]

    y[y < 0] = -y[y < 0]
    y[y > 100] = 200 - y[y > 100]

    # spremanje gustoće
    if k in save_steps:


        #podjela prostora na ćelije te broji cestice u svakoj celiji
        P, _, _ = np.histogram2d(  #broj cestica u celiji i,j, to je diskretna gustoća
            x, y,  #polozaji svih cestica
            bins=[nx, ny],  #dijelimo prostor na mrezu celija
            range=[[0, 100], [0, 100]]  #fizicki prostor
        )

#ovaj racuna gustocu
        #P je broj cestica u celiji ovisi o broju cestuica i velicini celije nije fizicka velicina 
        rho = P / (Nw * dx * dy)  #fizičku gustoću potrebnu za usporedbu s difuzijskom jednadžbom.
        results[k] = rho


#grid koordinata za graf
x_c = np.linspace(0, 100, nx)
y_c = np.linspace(0, 100, ny)
X, Y = np.meshgrid(x_c, y_c)



fig = plt.figure(figsize=(10, 7))

#crtam 4 grafa razlicitih vremena
for i, k in enumerate(save_steps):

    #2x2 raspored grafiva
    ax = fig.add_subplot(2, 2, i+1, projection='3d')

    Z = results[k].T   

    ax.plot_surface(X, Y, Z, cmap="Reds")

    ax.set_title("k = " + str(k))
    ax.set_xlabel("x")
    ax.set_ylabel("y")

plt.tight_layout()
plt.show()
