#### Série de Fourier, équation de la chaleur

## Modules

from math import *
from cmath import *
import matplotlib.pyplot as plt
import numpy as np

## Coefficients de Fourier et fonctions tests

'''
Convention : pour f une fonction 1-périodique,

c_n(f) := int_0^1 f(t)exp(-2 i pi n t) dt

Version discrète (par la méthode des rectangles) :

c_n^k(f) = (1/k)*sum_{l=0}^{k-1} f(l/k)exp(- 2 i pi l n / K).
'''

def fourier(f,k,n) : # renvoie c_n(f) approché avec k subdivisions, avec k un entier pair. On calcule la somme pour n allant de -k/2 à k/2
    #S1 = 0j
    #S2 = 0j
    S = 0
    for l in range(0,k) :
        z = -2j*pi*n*l/k
        S += f(l/k)*exp(z)
        #z1 = 2j*pi*l*n/k
        #z2 = -2j*pi*l*n/k
        #S1 += f(-l/k)*exp(z1)
        #S2 += f(l/k)*exp(z2)
    #return (S1+S2+f(0))/k
    return S/k

def poly_trigo(x) :
    return cos(6*pi*x)+0.5*sin(2*pi*x)
    #return cos(pi*x)
poly_trigo = np.vectorize(poly_trigo)

def creneau(x) :
    if x < 1/3 or x > 2/3 :
        return 0
    else :
        return 1
creneau = np.vectorize(creneau)

## Décomposistion en fréquence

k = 100 # subdivisions pour calculer les intégrales
N = 30
L_n = np.arange(-N,N+1)

L_poly = [fourier(poly_trigo,k,n) for n in L_n]
L_creneau = [fourier(creneau,k,n) for n in L_n]

plt.subplot(211)
plt.stem(L_n,np.abs(L_poly),markerfmt=' ',basefmt=' ')
plt.xlabel('Fréquence $n$')
plt.ylabel('$c_n(f)$')
plt.title('Décomposition en fréquence d un polynôme trigonométrique, K=%s'%k)
plt.show()

plt.subplot(212)
plt.stem(L_n,np.abs(L_creneau),markerfmt=' ',basefmt=' ')
plt.xlabel('Fréquence $n$')
plt.ylabel('$c_n(f)$')
plt.title('Décomposition en fréquence d un créneau, K=%s'%k)
plt.show()

## Comparaison avec la routine FFT ??

k = 100 # subdivisions pour calculer les intégrales
N = 100
L_n = np.arange(-N,N+1)

L_poly = [fourier(poly_trigo,k,n) for n in L_n]
L_creneau = [fourier(creneau,k,n) for n in L_n]

plt.subplot(221)
plt.stem(L_n,np.abs(L_poly),markerfmt=' ',basefmt=' ')
plt.xlabel('Fréquence $n$')
plt.ylabel('$c_n(f)$')
plt.title('Décomposition en fréquence d un polynôme trigonométrique, K=%s'%k)
plt.show()

plt.subplot(222)
plt.stem(L_n,np.abs(np.fft.fft(poly_trigo(L_n))),markerfmt=' ',basefmt=' ')
plt.xlabel('Fréquence $n$')
plt.ylabel('$c_n(f)$')
plt.title('Décomposition en fréquence d un polynôme trigonométrique via FFT, K=%s'%k)
plt.show()

plt.subplot(223)
plt.stem(L_n,np.abs(L_creneau),markerfmt=' ',basefmt=' ')
plt.xlabel('Fréquence $n$')
plt.ylabel('$c_n(f)$')
plt.title('Décomposition en fréquence d un créneau, K=%s'%k)
plt.show()

plt.subplot(224)
plt.stem(L_n,np.abs(np.fft.fft(L_creneau)),markerfmt=' ',basefmt=' ')
plt.xlabel('Fréquence $n$')
plt.ylabel('$c_n(f)$')
plt.title('Décomposition en fréquence d un créneau, via FFT K=%s'%k)
plt.show()

## Calcul de la somme partielle de Fourier : fonction polynômiale trigonométrique

def somme_fourier(f,x,N,k) : #renvoie la somme des c_n(f)(x)exp(...) pour n allant de -N/2 à N/2. On suppose que N est pair. k : nombre de discrétisation
    S1 = 0j
    S2 = 0J
    for n in range(1,int(N/2)) :
        c1 = fourier(f,k,n)
        z1 = 2j*pi*n*x
        S1 += c1*exp(z1)
    return 2*S1+fourier(f,k,0)


L_x = np.linspace(0,1,k)
k = 100

N = 2
plt.subplot(221)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

N = 4
plt.subplot(222)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

N = 6
plt.subplot(223)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

N = 10
plt.subplot(224)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

plt.show()

## Calcul de la somme partielle de Fourier : fonction créneau

def poly_trigo(x) : #flemme de remplacer
    #return creneau(x)
    return abs(x**2-1/2)*sqrt(x)*cos(4*pi*x)
poly_trigo = np.vectorize(poly_trigo)

L_x = np.linspace(0,1,k)
k = 100

N = 2
plt.subplot(221)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

N = 40
plt.subplot(222)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

N = 10
plt.subplot(223)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

N = 15
plt.subplot(224)
plt.plot(np.linspace(0,1,150),poly_trigo(np.linspace(0,1,150)).real,color='red',label='Fonction initiale')
S = [somme_fourier(poly_trigo,x,N,k).real for x in L_x]
plt.plot(L_x,S,'--',label='N=%s'%N)
plt.legend()

plt.show()

## Équation de la chaleur

























