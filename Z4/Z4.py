import numpy as np

# broj uzoraka
n1 = 5000000
n2 = 400000

# funkcija
def f(x):
    return np.exp(-x**2)  

#1.nacin
x1 = np.random.rand(n1)  #uzmes puno slucajnih br
Fn1 = np.mean(f(x1))     #prosjek = vrijednost integrala
sigma1 = np.sqrt(np.mean((f(x1) - Fn1)**2)) 



#sigma1 = np.std(f(x1))   #rasprsenost vrijednosti funk tj kolika je nesigurnost

#drugi nacin
A = np.exp(1)/(np.exp(1)-1)      #normalizacijska konstanta za p(x)=Ae^(-x)

u = np.random.rand(n2)    #slucajni br 
x2 = -np.log(1 - u*(1 - np.exp(-1)))   #pretvaram uniformne brojeve u željenu raspodjelu (jer su najvece vrijednosti blizu 0)

#povecavam i smanjujem utjecaj onih manje i vise biranih
weights = f(x2)/(A*np.exp(-x2)) #djeli sa p(x) - koliko cesto biram neki br

Fn2 = np.mean(weights)
# sigma2 = np.std(weights)
sigma2 = np.sqrt(np.mean((weights - Fn2)**2))  #mjeri nesigurnoist



with open("Z4.txt", "w") as f:
    f.write("p(x)=1:\n")
    f.write(f"Fn = {Fn1}\n")
    f.write(f"sigma = {sigma1}\n")
    f.write(f"sigma/sqrt(n) = {sigma1/np.sqrt(n1)}\n\n")

    f.write("p(x)=Ae^{-x}:\n")
    f.write(f"Fn = {Fn2}\n")
    f.write(f"sigma = {sigma2}\n")           #rasprsenost podataka
    f.write(f"sigma/sqrt(n) = {sigma2/np.sqrt(n2)}\n")    # nesigurnost srednje vrijednosti

