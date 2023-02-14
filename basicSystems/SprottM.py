#!/usr/bin/env python3

# Sprott-M

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D



def f(state, t):
    x, y, z = state  # Unpack the state vector
    
    xdot = -z
    ydot = -x**2 - y
    zdot = 1.7 + 1.7*x + y
    
    return xdot, ydot, zdot  # Derivatives


fig = plt.figure()
ax = fig.gca(projection="3d")


state0 = [0,0,0]
t = np.arange(0.0, 1000.0, 0.01)

states = odeint(f, state0, t)


ax.plot(states[1500:-1, 0], states[1500:-1, 1], states[1500:-1, 2],'b')
ax.plot(state0[0],state0[1],state0[2],'*')


ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')


plt.draw()
plt.show()