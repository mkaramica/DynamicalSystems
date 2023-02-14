#!/usr/bin/env python3

# Rounge-Kutta Solver for 2D space

import numpy as np
from pylab import *
import math

#Sample function:
    #   dx/dt=f(t,x,y)=x/t-ln(y)  ,  x0=0
    #   dy/dt=g(t,x,y)=(x+sint)y  ,  y0=1

#   x=tcost
#   y=exp(tsint)


def f(tFunc,xFunc,yFunc):
    return cos(tFunc)-math.log(yFunc,exp(1))

def g(tFunc,xFunc,yFunc):
    return (xFunc+sin(tFunc))*yFunc

def RK4(t0,x0,y0,dt,n):
    x_t=[x0]
    y_t=[y0]
    
    ti=t0
    
    for i in range(0,n-1):
        k1=f(ti,x_t[i],y_t[i])
        L1=g(ti,x_t[i],y_t[i])
        
        k2=f(ti+dt/2,x_t[i]+dt*k1/2,y_t[i]+dt*L1/2)
        L2=g(ti+dt/2,x_t[i]+dt*k1/2,y_t[i]+dt*L1/2)
        
        k3=f(ti+dt/2,x_t[i]+dt*k2/2,y_t[i]+dt*L2/2)
        L3=g(ti+dt/2,x_t[i]+dt*k2/2,y_t[i]+dt*L2/2)
        
        k4=f(ti+dt,x_t[i]+dt*k3,y_t[i]+dt*L3)
        L4=g(ti+dt,x_t[i]+dt*k3,y_t[i]+dt*L3)
        
        x_next=x_t[i]+1/6*dt*(k1+2*k2+2*k3+k4)
        y_next=y_t[i]+1/6*dt*(L1+2*L2+2*L3+L4)
        
        x_t.append(x_next)
        y_t.append(y_next)
        
        ti=ti+dt
    
    return x_t,y_t
    
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
x0=0    
y0=1    

tn=2*math.pi

n=2000

dt=(tn-t0)/n



t=linspace(t0, tn, num=n, endpoint=True)
x1=[]
y1=[]

[x2,y2]=RK4(t0,x0,y0,dt,n)
#x3=EulerMethod(t0,x0,dt,n)

for i in range(0,len(t)):
    x1.append(t[i]*cos(t[i]))
    y1.append(exp(t[i]*sin(t[i])))  


figure(1)
subplot(2,1,1)
plot(t,x1,'b',label = 'analytical')
plot(t,x2,'.r',label = 'Runge–Kutta')
title('x(t)')
grid(True)
legend(loc = 'best')

subplot(2,1,2)
plot(t,y1,'b',label = 'analytical')
plot(t,y2,'.r',label = 'Runge–Kutta')
title('y(t)')
grid(True)
legend(loc = 'best')

show()

figure(2)
plot(x1,y1,'b',label = 'analytical')
plot(x2,y2,'.r',label = 'Runge–Kutta')
grid(True)
legend(loc = 'best')

show()

#plot(t,x2,'.r',label = 'Runge–Kutta')
#plot(t,x3,'.g',label = 'Euler Method')

