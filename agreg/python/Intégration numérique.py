#### Intégration numérique

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Fonctions test

def f(x) :
    return exp(cos(x**3))

def p0(x) :
    return 1

def p1(x) :
    return x

def p2(x) :
    return x**2

def p3(x) :
    return x**3

def p4(x) :
    return x**4

alpha,beta = 0,1 # domaine d'intégration
L_I = [] # valeurs d'intégrales connues pour visualiser l'ordre des méthodes sur [0,1]
for i in range(0,5) :
    L_I.append(1/(i+1))

## Méthodes des rectangles à gauche et à droite

def rectangle_g(f,N,a,b) : #on prend ici une subdivision régulière
    S = 0
    dx = (b-a)/N # taille de la subdivision
    for i in range(N) :
        x_i = a + i*dx
        S += dx*f(x_i)
    return S

def rectangle_d(f,N,a,b) :
    S = 0
    dx = (b-a)/N
    for i in range(1,N+1) :
        x_i = a + i*dx
        S += dx*f(x_i)
    return S

Test_ordre = [abs(L_I[0]-rectangle_g(p0,30,0,1)),abs(L_I[1]-rectangle_g(p1,30,0,1)),abs(L_I[2]-rectangle_g(p2,30,0,1)),abs(L_I[3]-rectangle_g(p3,30,0,1)),abs(L_I[4]-rectangle_g(p4,30,0,1))]

ordre_rectangle = 0
for i in range(len(Test_ordre)-1) :
    if Test_ordre[i+1] <= 10**-5 :
        ordre_rectangle += 1

## Méthodes du point milieu

def point_milieu(f,N,a,b) :
    S = 0
    dx = (b-a)/N
    for i in range(N) :
        x_0 = a + i*dx
        x_1 = a + (i+1)*dx
        S += dx*f((x_0+x_1)/2)
    return S

Test_ordre = [abs(L_I[0]-point_milieu(p0,30,0,1)),abs(L_I[1]-point_milieu(p1,30,0,1)),abs(L_I[2]-point_milieu(p2,30,0,1)),abs(L_I[3]-point_milieu(p3,30,0,1)),abs(L_I[4]-point_milieu(p4,30,0,1))]

ordre_point_milieu = 0
for i in range(len(Test_ordre)-1) :
    if Test_ordre[i+1] <= 10**-5 :
        ordre_point_milieu += 1

## Méthode des trapèzes

def trapezes(f,N,a,b) :
    S = 0
    dx = (b-a)/N
    for i in range(N) :
        x_0 = a + i*dx
        x_1 = a + (i+1)*dx
        S += dx*(f(x_0)+f(x_1))/2
    return S

Test_ordre = [abs(L_I[0]-trapezes(p0,30,0,1)),abs(L_I[1]-trapezes(p1,30,0,1)),abs(L_I[2]-trapezes(p2,30,0,1)),abs(L_I[3]-trapezes(p3,30,0,1)),abs(L_I[4]-trapezes(p4,30,0,1))]

ordre_trapezes = 0
for i in range(len(Test_ordre)-1) :
    if Test_ordre[i+1] <= 10**-5 :
        ordre_trapezes += 1

## Méthode de Simpson

def simpson(f,N,a,b) :
    S = 0
    dx = (b-a)/N
    for i in range(N) :
        x_0 = a + i*dx
        x_1 = x_0 + dx/2
        x_2 = x_0 + dx
        S += ((1/6)*f(x_0) + (2/3)*f(x_1) + (1/6)*f(x_2))*dx
    return S

Test_ordre = [abs(L_I[0]-simpson(p0,30,0,1)),abs(L_I[1]-simpson(p1,30,0,1)),abs(L_I[2]-simpson(p2,30,0,1)),abs(L_I[3]-simpson(p3,30,0,1)),abs(L_I[4]-simpson(p4,30,0,1))]

ordre_simpson = 0
for i in range(len(Test_ordre)-1) :
    if Test_ordre[i+1] <= 10**-8 :
        ordre_simpson += 1

## Comparaison des vitesses de convergence

Iexp = np.exp(1) - np.exp(0) # intégrale de exp(x)dx sur [0,1]
def g(x) :
    return cos(20*x)
Ig = np.sin(20)/20

L_N = [N for N in np.arange(10,1000,100)]
I_rectangle_g = []
I_point_milieu = []
I_trapezes = []
I_simpson = []

'''
for N in L_N :
    I_rectangle_g.append(rectangle_g(exp,N,0,1))
    I_point_milieu.append(point_milieu(exp,N,0,1))
    I_trapezes.append(trapezes(exp,N,0,1))
    I_simpson.append(simpson(exp,N,0,1))

Erreur_rectangle_g = np.abs(Iexp - I_rectangle_g)
Erreur_point_milieu = np.abs(Iexp - I_point_milieu)
Erreur_trapezes = np.abs(Iexp - I_trapezes)
Erreur_simpson =  np.abs(Iexp - I_simpson)'''

for N in L_N :
    I_rectangle_g.append(rectangle_g(g,N,0,1))
    I_point_milieu.append(point_milieu(g,N,0,1))
    I_trapezes.append(trapezes(g,N,0,1))
    I_simpson.append(simpson(g,N,0,1))

Erreur_rectangle_g = np.abs(Ig - I_rectangle_g)
Erreur_point_milieu = np.abs(Ig - I_point_milieu)
Erreur_trapezes = np.abs(Ig - I_trapezes)
Erreur_simpson =  np.abs(Ig - I_simpson)

a = np.polyfit(np.log(L_N),np.log(Erreur_rectangle_g),1)[0]
b = np.polyfit(np.log(L_N),np.log(Erreur_point_milieu),1)[0]
c = np.polyfit(np.log(L_N),np.log(Erreur_trapezes),1)[0]
d = np.polyfit(np.log(L_N),np.log(Erreur_simpson),1)[0]

plt.plot(np.log(L_N),np.log(Erreur_rectangle_g),'-x',label='Rectangle à gauche, coeff =%f'%a)
plt.plot(np.log(L_N),np.log(Erreur_point_milieu),'-x',label='Point Milieu, coeff =%f'%b)
plt.plot(np.log(L_N),np.log(Erreur_trapezes),'-x',label='Trapezes, coeff =%f'%c)
plt.plot(np.log(L_N),np.log(Erreur_simpson),'-x',label='Simpson, coeff =%f'%d)

plt.legend()
plt.title("Courbes de convergence")
plt.xlabel('log(N)')
plt.ylabel('log(Erreur)')
plt.grid()
plt.show()


