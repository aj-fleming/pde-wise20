# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 14:11:48 2020

@author: AJ
"""

from meshops import MeshOperations
from prog2routines import setBoundaryValues, setEquationsSystemQuadElements

import matplotlib.pyplot as plt
import numpy as np

mesh = MeshOperations("mesh/quadrilateralmesh.msh")

basis = [
       lambda xi, eta: 1-xi-eta+xi*eta,
       lambda xi, eta: xi-xi*eta,
       lambda xi, eta: xi*eta,
       lambda xi, eta: eta-xi*eta
       ]
gradBasis = [
    lambda xi,eta: np.array([[-1+eta],[-1+xi]]),
    lambda xi,eta: np.array([[1-eta],[-xi]]),
    lambda xi,eta: np.array([[eta],[xi]]),
    lambda xi,eta: np.array([[-eta],[1-xi]])
    ]

def f(x,y):
    return 1

a1,b1 = setEquationsSystemQuadElements(mesh,basis,gradBasis,f)
a2,b2 = setBoundaryValues(mesh, a1.copy(), b1.copy(),1)

u = np.linalg.solve(a2,b2)

fig = plt.figure()
ax = fig.gca(projection='3d')

nodePos = mesh.getListOfNodePositions()

ax.plot_trisurf(nodePos[:,0].flatten(),nodePos[:,1].flatten(),u.flatten(),cmap='viridis')
plt.show()