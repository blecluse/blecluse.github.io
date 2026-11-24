#### Équation de transport 1D

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Différence finies à gauche

def DFG(J,T,u0,a,dt) :
    dx = 1/J # pas d'espace, on résout sur [0,1]
    N = int(T/dt) # nombre de points temps
    A = np.diag(np.ones(J))-np.diag(np.ones(J-1),k=-1)
    U = [u0(dx*j) for j in range(J)]
    for n in range(1,N) :
        U[0]=U[-1]
        U = U - a*dt/dx*np.dot(A,U)
    lam = a*dt/dx
    return U,lam

def DFD(J,T,u0,a,dt) :
    dx = 1/J # pas d'espace, on résout sur [0,1]
    N = int(T/dt) # nombre de points temps
    A = -np.eye(J) + np.eye(J,k=1)
    U = [u0(dx*j) for j in range(J)]
    for n in range(1,N) :
        U[0]=U[-1]
        U = U - a*dt/dx*np.dot(A,U)
    lam = a*dt/dx
    return U,lam

## Donnée initiale sinusoïdale

a = 1
J = 100 # nombre de points espace
T = 5 # temps de résolution
dt = 1/100 # pas de temps

L_x = np.linspace(0,1,J)

def u0(x) :
    return np.cos(6*pi*x)

sol_G = DFG(J,T,u0,a,dt)[0]
sol_D = DFD(J,T,u0,a,dt)[0]
lam = DFG(J,T,u0,a,dt)[1]

plt.plot(L_x,u0(L_x),label='Conditions initiales')
plt.plot(L_x,sol_G,label='Solutions approchés à t=%s, par DFG'%T)
#plt.plot(L_x,sol_D,label='Solutions approchés à t=%s, par DFD'%T)
plt.legend()
plt.show()

''' DFD fait n'imp '''

## Donnée initiale créneau

a = 0.9
J = 100 # nombre de points espace
T = 1 # temps de résolution
dt = 1/100 # pas de temps

L_x = np.linspace(0,1,J)

def u0(x) :
    if x <= 1/3 or x >= 2/3 :
        return 0
    else :
        return 1
u0 = np.vectorize(u0)

sol_G = DFG(J,T,u0,a,dt)[0]
lam = DFG(J,T,u0,a,dt)[1]

plt.plot(L_x,u0(L_x),label='Conditions initiales')
plt.plot(L_x,sol_G,label='Solutions approchés à t=%s, par DFG'%T)
plt.legend()
plt.show()

## Commentaires

'''

Différence finie à gauche : stable.
Différence finie à droite : instable.

On a posé lambda = a*dt/dx.
Si lambda = 1 : translation parfaite.
Si lambda > 1 : n'imp.
Si lambda < 0 : n'imp.
Si 0 < lambda < 1 : résultats convenables, on gagne de la régularité grâce à la diffusion numérique.

'''




