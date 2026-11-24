#### Résolution systèmes linéaires

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Matrices et solutions test

n = 5 # taille des matrices

A = 2*np.eye(n) - np.eye(n,k=-1) - np.eye(n,k=1) # matrice du Laplacien en dimension 1
B = 4*np.eye(n) + np.eye(n,k=-1) + np.eye(n,k=1) # matrice des fonctions splines

x = np.ones(n) # On choisit b, dans le système Ax=b, tel que x=(1, ... ,1)^T soit solution
bA = A.dot(x)
bB = B.dot(x)

## Programmation des méthodes

def decomposition(A) : # A = D + E + F
    n = np.shape(A)[0]
    D = np.diag(np.diag(A))
    F = np.zeros((n,n))
    for i in range(1,n) :
        F += np.diag(np.diag(A,k=i),k=i)
    E = A - D - F
    return D,E,F

def resolution(M,N,b,x0,eps,maxiter) :
    A = M - N
    i = 0 # Nombre d'itérations
    erreur = np.linalg.norm(A.dot(x0)-b) # Erreur
    x = x0
    L = [x] # Liste des valeurs
    while erreur > eps and i < maxiter :
        x = np.linalg.inv(M).dot(N).dot(x) + np.linalg.inv(M).dot(b)
        erreur = np.linalg.norm(A.dot(x)-b)
        i += 1
        L.append(x)
    return x,i,L

def Jacobi(A,b,x0,maxiter,eps) : # CV : A à diagonale strictement dominante
    D = decomposition(A)[0]
    M = D
    N = M - A
    return resolution(M,N,b,x0,eps,maxiter)

def Gauss_Seidel(A,b,x0,maxiter,eps) : # CV : A à diagonale strictement dominante ou A symétrique définie positive
    D,E,F = decomposition(A)
    M = E + D
    N = M - A
    return resolution(M,N,b,x0,eps,maxiter)

def relaxation(A,b,x0,maxiter,eps,omega) : # CV : pour A sysmétrique définie positive, CV ssi 0 < omega < 2
    D,E,F = decomposition(A)
    M = (1/omega)*D + E
    N = M - A
    return resolution(M,N,b,x0,eps,maxiter)

## Test des méthodes

x0 = np.random.randint(-10,10,n)
maxiter = 10**4
eps = 10**-5

# Jacobi

SolA1 = Jacobi(A,bA,x0,maxiter,eps)[0]
SolB1 = Jacobi(B,bB,x0,maxiter,eps)[0]

# Gauss-Seidel

SolA2 = Gauss_Seidel(A,bA,x0,maxiter,eps)[0]
SolB2 = Gauss_Seidel(B,bB,x0,maxiter,eps)[0]

# relaxation

w = 1
SolA3 = relaxation(A,bA,x0,maxiter,eps,w)[0]
SolB3 = relaxation(B,bB,x0,maxiter,eps,w)[0]

## Vitesse de convergence Jacobi

L_erreurJacobi_laplacien = [np.linalg.norm(x-xk) for xk in Jacobi(A,bA,x0,maxiter,eps)[2]]
L_erreurJacobi_splines = [np.linalg.norm(x-xk) for xk in Jacobi(B,bB,x0,maxiter,eps)[2]]
L_X_laplacien = range(len(Jacobi(A,bA,x0,maxiter,eps)[2]))
L_X_splines = range(len(Jacobi(B,bB,x0,maxiter,eps)[2]))

plt.plot(L_X_laplacien,np.log(L_erreurJacobi_laplacien),label='Matrice du Laplacien en dim 1, coeff')
plt.plot(L_X_splines,np.log(L_erreurJacobi_splines),label='Matrice des fonctions splines')
plt.legend()
plt.xlabel('Itérations')
plt.ylabel('log(Erreur)')
plt.title('Vitesse de convergence pour la méthode de Jacobi')

tauxAJacobi = exp(np.polyfit(L_X_laplacien,np.log(L_erreurJacobi_laplacien),1)[0])
tauxBJacobi = exp(np.polyfit(L_X_splines,np.log(L_erreurJacobi_splines),1)[0])
M = decomposition(A)[0]
N = M-A
RS_mat_iteration_A = max(np.linalg.eigvals(np.linalg.inv(M).dot(N)))
M = decomposition(B)[0]
N = M-B
RS_mat_iteration_B = max(np.linalg.eigvals(np.linalg.inv(M).dot(N)))

print('Taux de convergence estimé pour la matrice A et la méthode de Jacobi : ',tauxAJacobi)
print('Rayon spectral de la matrice d itération pour la matrice A : ',RS_mat_iteration_A)
print('Taux de convergence estimé pour la matrice B et la méthode de Jacobi : ',tauxBJacobi)
print('Rayon spectral de la matrice d itération pour la matrice B : ',RS_mat_iteration_B)
plt.show()

## Gauss-Seidel (idem)

## Relaxation (pour la matrice A)

L_omega = np.arange(0.01,2,0.02)
L_RS = []
D,E,F = decomposition(A)
for omega in L_omega :
    M = (1/omega)*D + E
    N = M - A
    RS = np.linalg.eigvals(np.linalg.inv(M).dot(N))
    L_RS.append(max(np.abs(RS)))

r = max(np.abs(np.linalg.eigvals(-np.linalg.inv(D).dot(E+F))))
omega0 = 2/(1+sqrt(1-r**2))


def rho(omega) :
    if 0 < omega < omega0:
        return (omega**2*r**2-2*(omega-1)+omega*r*sqrt(omega**2*r**2-4*omega+4))/2
    if omega0 <= omega < 2 :
        return omega-1
rho = np.vectorize(rho)

plt.close('all')
plt.plot(L_omega,L_RS,label='rho(omega) numériquement')
plt.plot(L_omega,rho(L_omega),'--',label='rho(omega) théorique')
plt.legend()
plt.ylabel('rho(omega)')
plt.xlabel('omega')
plt.show()



