#### Résolutions d'équations non linéaires - Méthodes itératives

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Fonctions tests

def f(x):
    return x**3-4*x+1
f = np.vectorize(f)

def df(x) :
    return 3*x**2-4
df = np.vectorize(df)

## Méthode de Newton 1-D

''' Critère de CV théorique : f de classe C^2 sur un intervalle I et f' ne s'annule pas sur cet intervalle.'''
''' Vitesse de convergence : en 2^p. |x_p-a| <= 1/C * [C*(|x_0-a|)]^(2^p), avec C une constante.'''

def Newton(f,df,x0,N) : # (fonction, fonction dérivée, valeur initiale, nombre d'itérations)
    x = x0
    for i in range(N) :
        x = x - f(x)/df(x)
    return x

L_N = np.arange(0,10,1)
L_x = np.arange(0.5,4,0.01)
L_0 = [0 for x in L_x]
L_f = f(L_x)
plt.grid()
plt.plot(L_x,L_f,label='y=x^3-4x+1')
plt.plot(L_x,L_0,'--',color='black')
L_Newton = [Newton(f,df,4,N) for N in L_N]
plt.plot(L_Newton,f(L_Newton),'x',color='red',label='x0=4')
plt.legend()
plt.title('Convergence de la méthode de Newton')
plt.show()

## Méthode de la sécante 1-D

def secante(f,df,x0,x1,N) : # (fonction, fonction dérivée, x0 et x1 encadrent la solution, nombre d'itérations)
    for i in range(N) :
        a = (f(x1)-f(x0))/(x1-x0)
        x0 = x1
        x1 = x1 - f(x1)/a
    return x1

L_N = np.arange(0,10,1)
L_x = np.arange(0.5,4,0.01)
L_0 = [0 for x in L_x]
L_f = f(L_x)
plt.grid()
plt.plot(L_x,L_f,label='y=x^3-4x+1')
plt.plot(L_x,L_0,'--',color='black')
L_secante = [secante(f,df,1,3,N) for N in L_N]
plt.plot(L_secante,f(L_secante),'x',color='red',label='x0=4')
plt.legend()
plt.title('Convergence de la méthode de la sécante')
plt.show()

## Vitesse de convergence des méthodes 1-D

def g(x) :
    return x**3+x
g = np.vectorize(g)

def dg(x) :
    return 3*x**2+1
dg = np.vectorize(dg)

Sol = 0 # unique solution de l'équation g(x)=0.

L_N = np.arange(1,8,1)
L_Newton = [Newton(g,dg,3,N) for N in L_N]
L_erreurN = [abs(x - Sol) for x in L_Newton] # Ici Sol=0 mais c'est pour le principe
L_log_erreurN = np.log(L_erreurN)
L_log_N = np.log(L_N)
plt.grid()
plt.plot(L_N,np.log(np.abs(L_log_erreurN)),label='Méthode de Newton')
L_secante = [secante(g,dg,-3,2,N) for N in L_N]
L_erreurS = [abs(x - Sol) for x in L_secante] # Ici Sol=0 mais c'est pour le principe
L_log_erreurS = np.log(L_erreurS)
plt.plot(L_N,np.log(np.abs(L_log_erreurS)),label='Méthode de la sécante')
plt.legend()
plt.ylabel('log(log(Erreur))')
plt.xlabel('Itérations')
plt.title('Comparaison des méthodes')
plt.show()

## Méthode de Newton-Raphson multi-D

'''Critère de CV théorique : f de classe C^2 et df inversible.'''

def F(X) :
    x,y = X[0],X[1]
    return np.array([x**2+x*y-2*y**2-4,x*exp(x)+y*exp(y)])

def dF(X) : # matrice 2x2
    x,y = X[0],X[1]
    return np.array([[2*x+y,x-4*y],[exp(x)*(1+x),exp(y)*(1+y)]])

def Newton_Raphson(F,dF,X0,N) :
    X = X0
    for i in range(N) :
        X = X - np.linalg.inv(dF(X)).dot(F(X))
    return X

L_N = np.arange(1,1000,10)
L_Newton_Raphson = [Newton_Raphson(F,dF,[2,2],N) for N in L_N]

#On peut de même illustrer la vitesse de convergence.'''





