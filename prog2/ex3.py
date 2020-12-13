# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:04:33 2020

@author: AJ
"""

import numpy as np
import matplotlib.pyplot as plt

from meshops import MeshOperations
from prog2routines import setEquationsSystem, setBoundaryValues
mesh = MeshOperations("mesh/unitSquareMany.msh")

basis2ndOrder = [
    lambda xi, eta: 1-3*(xi+eta)+2*(xi**2+eta**2)+4*xi*eta,
    lambda xi, eta:  xi*(2*xi-1),
    lambda xi, eta:  eta*(2*eta-1),
    lambda xi, eta:  4*xi*(1-xi-eta),
    lambda xi, eta:  4*xi*eta,
    lambda xi, eta:  4*eta*(1-xi-eta)
    
    ]
gradBasis2ndOrder = [
    lambda xi,eta: np.array([[-3+4*xi+4*eta],[-3+4*xi+4*eta]]),
    lambda xi,eta: np.array([[4*xi-1],[0]]),
    lambda xi,eta: np.array([[0],[4*eta-1]]),
    lambda xi,eta: np.array([[4-8*xi-4*eta],[-4*xi]]),
    lambda xi,eta: np.array([[4*eta],[4*xi]]),
    lambda xi,eta: np.array([[-4*eta],[4-8*eta-4*xi]])
    ]


a,b = setEquationsSystem(mesh,basis2ndOrder, gradBasis2ndOrder, 2, lambda x,y: np.sin(16*(x-0.5)*(y-0.5)*np.pi))
a, b = setBoundaryValues(mesh, a, b, 2)
u = np.linalg.solve(a,b)

fig = plt.figure()
ax = fig.gca(projection='3d')
mesh.plot(ax,u)
plt.show()