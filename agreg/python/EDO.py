#### EDO

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Euler explicite

def Euler_explicite(Y0,f,N,T) :
    dt = T/N # pas de la subdivision
    t = 0
    Y = Y0
    S = [Y]
    for i in range(N):
        Y = Y + dt*f(Y,t)
        t += dt
        S.append(Y)
    return S

## Euler implicite

def Newton(f,df,x0,eps) :
    x = x0
    while abs(f(x)) > eps :
        x = x - f(x)/df(x)
    return x

def Euler_implicite(Y0,f,df,N,T) : # df la dérivée de f par rapoprt à la variable d'état
    dt = T/N # pas de la subdivision
    t = 0
    Y = Y0
    S = [Y]
    for i in range(N):
        def g(x) :
            return x-Y-dt*f(t,x)
        def dg(x) :
            return 1-dt*df(t,x)
        Y = Newton(g,dg,Y,10**-10)
        t += dt
        S.append(Y)
    return S

## Test Euler explicite

def f(y,t) :
    return -4*t*y

def solution(t) :
    return np.exp(-2*t**2)

T = 2
L_x = np.linspace(0,T,100)

for N in [20] :
    S = Euler_explicite(1,f,N,T)
    L_x_N = np.linspace(0,T,N+1)
    plt.plot(L_x_N,S,'-x',label='Solution approchée, N=%s'%N)
    plt.legend()
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=-4t$, $y(0)=1$ avec Euler explicite')
plt.show()

## Test Euler implicite

def f(y,t) :
    return -4*t*y

def df(y,t) :
    return -4*t

def solution(t) :
    return np.exp(-2*t**2)

T = 2
L_x = np.linspace(0,T,100)

for N in [20] :
    S = Euler_implicite(1,f,df,N,T)
    L_x_N = np.linspace(0,T,N+1)
    plt.plot(L_x_N,S,'-x',label='Solution approchée, N=%s'%N)
    plt.legend()
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=-4t$, $y(0)=1$ avec Euler implicite')
plt.show()

## Comparaison implicite et explicite sur un problème simple

def f(y,t) :
    return -4*t*y

def df(y,t) :
    return -4*t

def solution(t) :
    return np.exp(-2*t**2)

T = 2
L_x = np.linspace(0,T,100)

for N in [20] :
    S_i = Euler_implicite(1,f,df,N,T)
    S_e = Euler_explicite(1,f,N,T)
    L_x_N = np.linspace(0,T,N+1)
    plt.plot(L_x_N,S_i,'--',label='Euler implicite, N=%s'%N)
    plt.plot(L_x_N,S_e,'--',label='Euler explicite, N=%s'%N)
    plt.legend()
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=-4t$, $y(0)=1$ avec Euler implicite')
plt.xlabel('Temps')
plt.show()

## Vitesse de convergence

L_N = [100*N for N in range(2,20)]
Erreur_e = []
Erreur_i = []

def f(y,t) :
    return -4*t*y

def df(y,t) :
    return -4*t

def solution(t) :
    return np.exp(-2*t**2)

T = 2
L_x = np.linspace(0,T,100)

for N in L_N :
    L_x_N = np.linspace(0,T,N+1)
    S_i = Euler_implicite(1,f,df,N,T)
    S_e = Euler_explicite(1,f,N,T)
    Solution_exacte = solution(L_x_N)
    Erreur_e.append(max(np.abs(Solution_exacte-S_e)))
    Erreur_i.append(max(np.abs(Solution_exacte-S_i)))
e = np.polyfit(np.log(L_N),np.log(Erreur_e),1)[0]
i = np.polyfit(np.log(L_N),np.log(Erreur_i),1)[0]
plt.plot(np.log(L_N),np.log(Erreur_e),label='Euler explicite, coefficient directeur = %s'%e)
plt.plot(np.log(L_N),np.log(Erreur_i),label='Euler implicite, coefficient directeur = %s'%i)
plt.title('Comparaison des vitesses de convergence')
plt.legend()
plt.show()


## Problème raide, comparaison Euler explicite et Euler implicite ????

def f(y,t) :
    return a*y

def df(y,t) :
    return a

def solution(t) :
    return np.exp(a*t)

T = 3
L_x = np.linspace(0,T,100)
N = 400

a = 1
S_i = Euler_implicite(1,f,df,N,T)
S_e = Euler_explicite(1,f,N,T)
L_x_N = np.linspace(0,T,N+1)
plt.subplot(221)
plt.plot(L_x_N,S_i,label='Euler implicite, N=%s'%N)
plt.plot(L_x_N,S_e,label='Euler explicite, N=%s'%N)
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=ay$, $y(0)=1$, a=%s'%a)
plt.xlabel('Temps')

a = 2
S_i = Euler_implicite(1,f,df,N,T)
S_e = Euler_explicite(1,f,N,T)
L_x_N = np.linspace(0,T,N+1)
plt.subplot(222)
plt.plot(L_x_N,S_i,label='Euler implicite, N=%s'%N)
plt.plot(L_x_N,S_e,label='Euler explicite, N=%s'%N)
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=ay$, $y(0)=1$, a=%s'%a)
plt.xlabel('Temps')

a = 4
S_i = Euler_implicite(1,f,df,N,T)
S_e = Euler_explicite(1,f,N,T)
L_x_N = np.linspace(0,T,N+1)
plt.subplot(223)
plt.plot(L_x_N,S_i,label='Euler implicite, N=%s'%N)
plt.plot(L_x_N,S_e,label='Euler explicite, N=%s'%N)
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=ay$, $y(0)=1$, a=%s'%a)
plt.xlabel('Temps')

a = 10
S_i = Euler_implicite(1,f,df,N,T)
S_e = Euler_explicite(1,f,N,T)
L_x_N = np.linspace(0,T,N+1)
plt.subplot(224)
plt.plot(L_x_N,S_i,label='Euler implicite, N=%s'%N)
plt.plot(L_x_N,S_e,label='Euler explicite, N=%s'%N)
plt.plot(L_x,solution(L_x),label='Solution exacte')
plt.legend()
plt.title('Résolution $dy/dt=ay$, $y(0)=1$, a=%s'%a)
plt.xlabel('Temps')









plt.show()