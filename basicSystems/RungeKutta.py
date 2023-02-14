#!/usr/bin/env python3

# Rounge-Kutta Solver

import numpy as np
from pylab import *
import math

#Sample function dx/dt=x*sint, x0=1
#   x=exp(1-cost)

def f(tFunc,xFunc):
    return xFunc*sin(tFunc)

def RK4(t0,x0,dt,n):
    x_t=[x0]
    ti=t0
    
    for i in range(0,n-1):
        k1=f(ti,x_t[i])
        k2=f(ti+dt/2,x_t[i]+dt*k1/2)
        k3=f(ti+dt/2,x_t[i]+dt*k2/2)
        k4=f(ti+dt,x_t[i]+dt*k3)
        
        x_next=x_t[i]+1/6*dt*(k1+2*k2+2*k3+k4)
        x_t.append(x_next)
        ti=ti+dt
    
    
    return x_t
    
def EulerMethod(t0,x0,dt,n): 
    x_t=[x0]
    ti=t0
    
    for i in range(0,n-1):
        k1=f(ti,x_t[i])
        k2=f(ti+dt/2,x_t[i]+dt*k1/2)
        k3=f(ti+dt/2,x_t[i]+dt*k2/2)
        k4=f(ti+dt,x_t[i]+dt*k3)
        
        x_next=x_t[i]+dt*f(ti,x_t[i])
        x_t.append(x_next)
        ti=ti+dt
    return x_t
    
    

t0=0
x0=1    

tn=4*math.pi

n=200

dt=(tn-t0)/n



t=linspace(t0, tn, num=n, endpoint=True)
x1=[]

x2=RK4(t0,x0,dt,n)
x3=EulerMethod(t0,x0,dt,n)

for i in range(0,len(t)):
    x1.append(exp(1-cos(t[i])))

plot(t,x1,'b',label = 'analytical')
plot(t,x2,'.r',label = 'Runge–Kutta')
plot(t,x3,'.g',label = 'Euler Method')
title('Solution for 1st order differential equation: dx/dt=x*sint with x0=1')
grid(True)
legend(loc = 'best')
show()