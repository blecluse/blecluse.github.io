#### Équations des ondes

## Modules

from math import *
import matplotlib.pyplot as plt
import numpy as np

## Schéma différences finies

def DF(w,v,T,dx,dt,c) : # on résout pour x \in  [0,1]
    N = int(T/dt)
    J = int(1/dx)
    a = (c*dt/dx)**2
    A = -2*np.eye(J)+np.eye(J,k=1)+np.eye(J,k=-1)
    U0 = np.array([w(j*dx) for j in range(J)])
    U1 = U0+dt*np.array([v(j*dx) for j in range(J)])
    U0[0]=U0[-1]=0
    U1[0]=U1[-1]=0
    for n in range(2,N) :
        U = 2*U1-U0+a*np.dot(A,U1)
        U[0]=U[-1]=0
        U0 = U1
        U1 = U
    return U

## Tests

''' Conditions CFL : a <= 1'''

def u0(x) :
    return sin(2*pi*x)

def u1(x) :
    return 0

T = 5
dx = 0.01
dt = 0.9*dx
c = 1
a = (c*dt/dx)**2

def solution(x,t) :
    return ( sin(2*pi*(x+c*t))+sin(2*pi*(x-c*t)) + (cos(2*pi*(x-c*t))-cos(2*pi*(x+c*t)))/(2*c*pi) )/2
solution = np.vectorize(solution)

L_x = np.linspace(0,1,int(1/dx))
Sol = DF(u0,u1,T,dx,dt,c)
Sol_exacte = [solution(x,T) for x in np.linspace(0,1,100)]
plt.plot(np.linspace(0,1,100),Sol_exacte,label='Solution exacte')
plt.plot(L_x,Sol,label='Solution approchée, T=%s'%T)
plt.legend()
plt.title('Équations des ondes')
plt.show()