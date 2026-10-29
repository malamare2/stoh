import numpy as np
import matplotlib.pyplot as plt


No_s = 10000
No_t = 0
S_polut = 2.5
T_polut = 12.5

#konstanta raspada
lambda_S = np.log(2)/S_polut   # veza poluraspada i konstante raspada, koliko brzo se raspada
lambda_T = np.log(2)/T_polut

dt = 0.1    #gledama rospad svaki 0.1 min tj svakih 6sek koliko se raspalo
t_max = 100 # gledamo raspad od 0.1min do 100min
#stavaramo niz vremena 
times = np.arange(0, t_max, dt)  #ovo je lista svih vremenskih perioda koje promatramo od 0 do 100min

N_s = np.zeros(len(times))  # trazi s u trenutku times(i) - prazne liste za spremanje podataka
N_t = np.zeros(len(times))

N_s[0] = No_s  #početni prvi clan
N_t[0] = No_t

#svaka jezgra ima vjerojatnost da se raspadne u koraku
for i in range(1, len(times)):   #idemo kroz svaki ttrenutak u vremenu
    Sb = int(N_s[i-1]) #uzima trenutno stanje sustava iz proslog koraka i pretvaran u cijeli broj 
    Te = int(N_t[i-1])

    # raspad Sb → Te
    decay_Sb = np.sum(np.random.rand(Sb) < lambda_S * dt)  #sum gleda broj koliko ih se raspalo, uzima random br od sb tj jezgre i gleda jesu se raspale

    # raspad Te → dalje (gubitak)
    decay_Te = np.sum(np.random.rand(Te) < lambda_T * dt)

    N_s[i] = Sb - decay_Sb     #gube se jezgre koje su se raspale pa ih moramo oduzet
    N_t[i] = Te + decay_Sb - decay_Te   #prvo nastaju iz sb pa se raspadaju te ih oduzimamo



t = times

#eksponencijalni raspad
N_Sb_analytical = No_s * np.exp(-lambda_S * t)  #pocinjem od 10000 broj pada eksponencijalno što je lambda veći to brže pada

#Batemanova jednadžba za ulančani raspad - koliko nastaje pa koliko se rspada
N_Te_analytical = (
    No_s * lambda_S / (lambda_T - lambda_S)      #samo “skalira” koliko brzo Sb proizvodi Te i koliko brzo Te nestaje
    * (np.exp(-lambda_S * t) - np.exp(-lambda_T * t))  #koliko Sb jos postoji i koliko se Te samo raspada
)


plt.figure()
plt.plot(t, N_s, label="Sb (MC)")
plt.plot(t, N_t, label="Te (MC)")

plt.plot(t, N_Sb_analytical, "--", label="Sb analitički")
plt.plot(t, N_Te_analytical, "--", label="Te analitički")

plt.xlabel("vrijeme")
plt.ylabel("broj jezgri")
plt.legend()
plt.grid()

plt.show()