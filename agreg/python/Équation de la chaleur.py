#### Équation de la chaleur

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Différences finies

dt = 0.1 # discrétisation en temps
dx = 0.1 # discrétisation en espace
L = pi # longueur du segment, valeur fixée
T = 10 # temps de résolution
J = pi/dx # nombre de point en espace
N = T/dt # nombre de point en temps
a = dt/dx**2 # coefficient CFL

def DF(u0,dt,dx,T) :
    N = int(T/dt)
    U0 = np.array([u0(j*dx) for j in range(J)])
    U = U0
    a = dt/dx**2
    A = (1-2*a)*np.eye(J)+a*(np.eye(J,k=1)+np.eye(J,k=-1))
    print(A)
    for n in range(1,N) :
        U[0]=U[-1]=0
        U = np.dot(A,U)
    a = dt/dx**2
    return U,a

def dirac(x) :
    if pi/2-0.01 < x < pi/2 + 0.01 :
        return 100
    else :
        return 0
dirac = np.vectorize(dirac)

def gauss(x) :
    return 100*exp((-(x-pi/2)**2)/0.1)
gauss = np.vectorize(gauss)


'''
Le schéma est stable ssi a <= 1/2. On distingue donc trois cas.
'''

## Schéma stable : a < 1/2, pas de temps très petit

dt = 0.00001
dx = 0.01
L = pi # longueur du segment, valeur fixée
T = 0.01 # temps de résolution
J = int(pi/dx) # nombre de point en espace
N = int(T/dt) # nombre de point en temps
a = dt/dx**2 # coefficient CFL

L_x = np.linspace(0,pi,J)
L_U = DF(gauss,dt,dx,T)[0]

plt.subplot(221)
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(222)
T = 0.01
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(223)
T = 0.02
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(224)
T = 0.1
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.show()

## Schéma instable : a > 1/2

dt = 0.01
dx = 0.01
L = pi # longueur du segment, valeur fixée
T = 0.1 # temps de résolution
J = int(pi/dx) # nombre de point en espace
N = int(T/dt) # nombre de point en temps
a = dt/dx**2 # coefficient CFL

def gauss(x) : # pas une gaussienne mais flemme de changer le nom partout
    return pi/2-abs(x-pi/2)

L_x = np.linspace(0,pi,J)
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Conditions initiales, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()
plt.title('Méthode des différences finies à t=%s'%T)
plt.show()

## Cas limite : a = 1/2

dt = 0.00001
dx = sqrt(2*dt)
L = pi # longueur du segment, valeur fixée
T = 0.1 # temps de résolution
J = int(pi/dx) # nombre de point en espace
N = int(T/dt) # nombre de point en temps
a = dt/dx**2 # coefficient CFL

L_x = np.linspace(0,pi,J)
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Conditions initiales, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()
plt.title('Méthode des différences finies à t=%s'%T)
plt.show()

## Autre condition initiale

def gauss(x) : # pas une gaussienne mais flemme de changer le nom partout
    return pi/2-abs(x-pi/2)

dt = 0.00001
dx = 0.01
L = pi # longueur du segment, valeur fixée
T = 0.01 # temps de résolution
J = int(pi/dx) # nombre de point en espace
N = int(T/dt) # nombre de point en temps
a = dt/dx**2 # coefficient CFL

L_x = np.linspace(0,pi,J)
L_U = DF(gauss,dt,dx,T)[0]

plt.subplot(221)
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(222)
T = 0.01
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(223)
T = 0.5
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(224)
T = 1
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.show()

## Autre condition initiale

def gauss(x) : # pas une gaussienne mais flemme de changer le nom partout
    if pi/2-0.1 < x < pi/2+0.1 :
        return 1
    if pi/4-0.1 < x < pi/4+0.1 :
        return 1
    if 3*pi/4-0.1 < x < 3*pi/4+0.1 :
        return 1
    else :
        return 0
gauss = np.vectorize(gauss)

dt = 0.00001
dx = 0.01
L = pi # longueur du segment, valeur fixée
T = 0.01 # temps de résolution
J = int(pi/dx) # nombre de point en espace
N = int(T/dt) # nombre de point en temps
a = dt/dx**2 # coefficient CFL

L_x = np.linspace(0,pi,J)
L_U = DF(gauss,dt,dx,T)[0]

plt.subplot(221)
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(222)
T = 0.001
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(223)
T = 0.005
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.subplot(224)
T = 0.01
L_U = DF(gauss,dt,dx,T)[0]
plt.plot(L_x,gauss(L_x),label='Condition initiale, t=0')
plt.plot(L_x,L_U,label='Solution à t=%s'%T)
plt.ylabel('Température')
plt.xlabel('Espace')
plt.legend()

plt.show()