#!/usr/bin/env python3

# biforcation diagram for logistic mapping

import numpy as np
from pylab import *


Amin=1
Amax=3.8
dA=0.001



Amat=[]
Xmat=[]

ARow=[]

Qcurve1=[]
Qcurve2=[]
Qcurve3=[]
Qcurve4=[]
Qcurve5=[]
Qcurve6=[]

def logisticFucntion(Ain,xn):
    return xn+Ain*sin(xn)
    
def fof_logistic(n,Ain,xn):
    y=logisticFucntion(Ain,xn)
    for i in range(1,n):
        y=logisticFucntion(Ain,y)
    return y
        
    

nsec=int((Amax-Amin)/dA)


for i in range(0,nsec):

    A=Amin +dA*i
    x0=0.5
    ARow.append(A)
    Qcurve1.append(fof_logistic(1,A,0.5))
    Qcurve2.append(fof_logistic(2,A,0.5))
    Qcurve3.append(fof_logistic(3,A,0.5))
    Qcurve4.append(fof_logistic(4,A,0.5))
    Qcurve5.append(fof_logistic(5,A,0.5))
    Qcurve6.append(fof_logistic(6,A,0.5))
    
    for j in range(0,1000):
        x1=logisticFucntion(A,x0)
        x0=x1
        if j > 950:
            Amat.append(A)
            Xmat.append(x1)
figure(1)
plot(Amat,Xmat,'.',label = 'biforcation')
#plot(ARow,Qcurve1,'b',label = 'Q-Curve#1')
#plot(ARow,Qcurve2,'r',label = 'Q-Curve#2')
#plot(ARow,Qcurve3,'g',label = 'Q-Curve#3')
#plot(ARow,Qcurve4,'k',label = 'Q-Curve#4')
#plot(ARow,Qcurve5,'m',label = 'Q-Curve#5')
#plot(ARow,Qcurve6,'y',label = 'Q-Curve#6')
legend(loc = 'lower left')
xlim(Amin,Amax)
xlabel('A')
ylabel('X_inf')
grid(True)
title('Biforcation Diagram')
show()

exit()

x=linspace(0, 1, num=1000, endpoint=True)
y=[]
for i in range(0,len(x)):
    y.append(fof_logistic(nFoFselected,A_selected,x[i]))
    

figure(2)
plot(x,x,'--r')
plot(x,y,'b')
xlim(0,1)
xlabel('x_n')
ylabel('x_n+1')
grid(True)
title('Logistic n_th Mapping for n='+str(nFoFselected)+', A='+str(A_selected))
show()

# Butterfly Effect

A=4
x0=0.2
y0=x0+1e-8

xMat=[x0]
yMat=[y0]
tMat=[0]

for i in range(1,100):
    x0=A*x0*(1-x0)
    y0=A*y0*(1-y0)
    xMat.append(x0)
    yMat.append(y0)
    tMat.append(i)

figure(3)
plot(tMat,xMat,'*--r', label = 'x0=0.2')
plot(tMat,yMat,'.--b', label = 'x0=0.2000002')
xlabel('n')
ylabel('x_n')
grid(True)
title('Butterfly Effect for Logistic Equation , A=4')
legend(loc = 'lower left')
show()
      

exit()
    
