# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from meshops import MeshOperations
from prog2routines import setEquationsSystem, setBoundaryValues

   
basis1stOrder = [
       lambda xi, eta: 1-xi-eta,
       lambda xi, eta: xi,
       lambda xi, eta: eta
       ]
gradBasis1stOrder = [
    lambda xi,eta: -1*np.ones((2,1)),
    lambda xi,eta: np.array([[1],[0]]),
    lambda xi,eta: np.array([[0],[1]])
    ]


def f1c(x,y):
    return np.sin(2*np.pi*x)

def f1a(x,y):
    return 1
  
    
mesh_A = MeshOperations("mesh/unitSquare1.msh")
mesh_B = MeshOperations("mesh/unitSquare2.msh")
order = 1

A1a, b1a = setEquationsSystem(mesh_A, basis1stOrder, gradBasis1stOrder, order, f1a)
A1b, b1b = setEquationsSystem(mesh_B, basis1stOrder, gradBasis1stOrder, order, f1a)
A1c, b1c = setEquationsSystem(mesh_B, basis1stOrder, gradBasis1stOrder, order, f1c)

A1a2, b1a2 = setBoundaryValues(mesh_A, A1a.copy(), b1a.copy(), order)
A1b2, b1b2 = setBoundaryValues(mesh_B, A1b.copy(), b1b.copy(), order)
A1c2, b1c2 = setBoundaryValues(mesh_B, A1c.copy(), b1c.copy(), order)

u1a = np.linalg.solve(A1a2, b1a2)
u1b = np.linalg.solve(A1b2, b1b2)
u1c = np.linalg.solve(A1c2, b1c2)

fig = plt.figure(1)
fig.suptitle("Solution to 1a: $-\Delta u = 1$ on unitSquare1")
mesh_A.plot(fig.gca(projection='3d'), u1a)
fig = plt.figure(2)
fig.suptitle("Solution to 1b: $-\Delta u = 1$ on unitSquare2")
mesh_B.plot(fig.gca(projection='3d'), u1b)
fig = plt.figure()
fig.suptitle("Solution to 1c: $-\Delta u = \sin(2*\pi*x)$ on unitSquare2")
mesh_B.plot(fig.gca(projection='3d'), u1c)
plt.show()
    




