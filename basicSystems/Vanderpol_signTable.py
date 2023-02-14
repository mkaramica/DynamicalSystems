#!/usr/bin/env python3



import numpy as np
from pylab import *
import math



mio=3

def f(tFunc,xFunc,yFunc):
    return yFunc

def g(tFunc,xFunc,yFunc):
    return mio*yFunc*(1-xFunc*xFunc)-xFunc
    
  
    

def RK4(t0,x0,y0,dt,n):
    x_t=[x0]
    y_t=[y0]
    
    ti=t0
    
    for i in range(0,n-1):
        k1=f(ti,x_t[i],y_t[i])
        L1=f(ti,x_t[i],y_t[i])
        
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
    
def EulerMethod(t0,x0,y0,dt,n): 
    x_t=[x0]
    y_t=[y0]
    
    ti=t0
    
    for i in range(0,n-1):      
        x_next=x_t[i]+dt*f(ti,x_t[i],y_t[i])
        y_next=y_t[i]+dt*g(ti,x_t[i],y_t[i])
        x_t.append(x_next)
        y_t.append(y_next)
        ti=ti+dt
    return x_t,y_t
    
    

t0=0
x0=2
y0=4

tn=100

n=5000

dt=(tn-t0)/n



t=linspace(t0, tn, num=n, endpoint=True)
x1=[]
y1=[]

[x2,y2]=RK4(t0,x0,y0,dt,n)
[x3,y3]=EulerMethod(t0,x0,y0,dt,n)

for i in range(0,len(t)):
    x1.append(t[i]*cos(t[i]))
    y1.append(exp(t[i]*sin(t[i])))  




#Extremum Points:

xplot=linspace(-3, 3, num=1005, endpoint=True)

yExtremum=[]
yAlternative=[]

for i in range(0,len(xplot)):
    yExtremum.append(xplot[i]/(mio*(1-xplot[i]**2)))


tr_j=[]
for i in range(len(t)):
    tr_j.append(mio*(1-x2[i]**2))



figure(1)
#plot(x1,y1,'b',label = 'analytical')
plot(x2,y2,'b',label = 'Runge–Kutta')
#plot(x3,y3,'r',label = 'Euler')

plot(xplot,yExtremum,'g', label='Extremum-Y')



grid(True)
legend(loc = 'best')
title('Vanderpole System for mio=3')
xlabel('X')
ylabel('Y')
xlim(-3,3)
ylim(-6,6)
show()

figure(2)
plot(t,tr_j)

show()



