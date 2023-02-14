#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import math

m=1
k=.5

theta=math.pi/6

def f(state, t):
    x, y = state  # Unpack the state vector
    
    xdot=y
    ydot=-k/m*x**3
    
    
    return xdot, ydot  # Derivatives

state0 = [1.0, 0]
t = np.arange(0.0, 200.0, 0.01)

states = odeint(f, state0, t)

fig = plt.figure()
ax = fig.gca()
ax.plot(states[:, 0], states[:, 1],'r')
ax.plot(state0[0],state0[1],'*r')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.grid(True)
plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.gca().set_aspect('equal', adjustable='box')
plt.draw()
plt.show()