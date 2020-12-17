# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:53:39 2020

@author: AJ
"""

from meshops import MeshOperations
from prog2routines import setBoundaryValues, setEquationsSystem
import numpy as np
import matplotlib.pyplot as plt

order = 1
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

def f1a(x,y):
    return 1

mesh64 = MeshOperations("mesh/mesh64.msh")
mesh256 = MeshOperations("mesh/mesh256.msh")
mesh1024 = MeshOperations("mesh/mesh1024.msh")

a64, b64 = setBoundaryValues(mesh64,
                             *setEquationsSystem(mesh64,basis1stOrder,gradBasis1stOrder,order,f1a),
                             order)
a256, b256 = setBoundaryValues(mesh256,
                             *setEquationsSystem(mesh256,basis1stOrder,gradBasis1stOrder,order,f1a),
                             order)
a1024, b1024 = setBoundaryValues(mesh1024,
                             *setEquationsSystem(mesh1024,basis1stOrder,gradBasis1stOrder,order,f1a),
                             order)
u64=np.linalg.solve(a64,b64)
u256 = np.linalg.solve(a256,b256)
u1024 = np.linalg.solve(a1024,b1024)

fig = plt.figure(1)
fig.suptitle("Solution to $\Delta u = 1$ on mesh64")
mesh64.plot(fig.gca(projection='3d'), u64)
fig = plt.figure(2)
fig.suptitle("Solution to $\Delta u = 1$ on mesh256")
mesh256.plot(fig.gca(projection='3d'), u256)
fig = plt.figure(3)
fig.suptitle("Solution to $\Delta u = 1$ on mesh1024")
mesh1024.plot(fig.gca(projection='3d'), u1024)

errs = [
            mesh64.calcL2ErrorPoisson(u64,order,basis1stOrder),
            mesh256.calcL2ErrorPoisson(u256,1,basis1stOrder),
            mesh1024.calcL2ErrorPoisson(u1024,1,basis1stOrder)
        ]
fig = plt.figure(4)
ax = fig.gca()
fig.suptitle("$L^2$ error in solutions to $\Delta u = 1$")
ax.set_xscale("log",basex=2)
ax.set_xlabel("Number of Elements")
ax.set_yscale("log")
ax.set_ylabel("Measured $L^2$ Error")
ax.plot([64,256,1024],errs)
ax.scatter([64,256,1024], errs, marker="+",color="red")

plt.show()