# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from meshops import MeshOperations, triangleIntegrationRule
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

basis2ndOrder = [
    lambda xi, eta: 1-3*(xi+eta)+2*(xi**2+eta**2)+4*xi*eta,
    lambda xi, eta:  xi*(2*xi-1),
    lambda xi, eta:  eta*(2*eta-1),
    lambda xi, eta:  4*xi*(1-xi-eta),
    lambda xi, eta:  4*eta*(1-xi-eta),
    lambda xi, eta:  4*xi*eta
    ]
gradBasis2ndOrder = [
    lambda xi,eta: np.array([[-3+4*xi+4*eta],[-3+4*xi+4*eta]]),
    lambda xi,eta: np.array([[4*xi-1],[0]]),
    lambda xi,eta: np.array([[0],[4*eta-1]]),
    lambda xi,eta: np.array([[4-8*xi-4*eta],[-4*xi]]),
    lambda xi,eta: np.array([[-4*eta],[4-8*eta-4*xi]]),
    lambda xi,eta: np.array([[4*eta],[4*xi]])
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
fig.suptitle("Solution to 1a: f(x,y) = 1 on unitSquare1")
mesh_A.plot(fig.gca(projection='3d'), u1a)
fig = plt.figure(2)
fig.suptitle("Solution to 1b: f(x,y) = 1 on unitSquare2")
mesh_B.plot(fig.gca(projection='3d'), u1b)
fig = plt.figure(3)
fig.suptitle("Solution to 1c: f(x,y) = sin(2*pi*x) on unitSquare2")
mesh_B.plot(fig.gca(projection='3d'), u1c)

    
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

fig = plt.figure(4)
fig.suptitle("Solution to $\Delta u = 1$ on mesh64")
mesh64.plot(fig.gca(projection='3d'), u64)
fig = plt.figure(5)
fig.suptitle("Solution to $\Delta u = 1$ on mesh256")
mesh256.plot(fig.gca(projection='3d'), u256)
fig = plt.figure(6)
fig.suptitle("Solution to $\Delta u = 1$ on mesh1024")
mesh1024.plot(fig.gca(projection='3d'), u1024)

errs = [
            mesh64.calcL2ErrorPoisson(u64,order,basis1stOrder),
            mesh256.calcL2ErrorPoisson(u256,1,basis1stOrder),
            mesh1024.calcL2ErrorPoisson(u1024,1,basis1stOrder)
        ]
plt.show()



