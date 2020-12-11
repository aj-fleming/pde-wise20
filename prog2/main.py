# -*- coding: utf-8 -*-

from meshops import MeshOperations, triangleIntegrationRule
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def setEquationsSystem(meshdata, phi, delphi, order, forcingfunc):
    global_A = np.zeros((meshdata.getNodeCount(),meshdata.getNodeCount()))
    global_b = np.zeros((meshdata.getNodeCount(),1))
    
    [quadp, quadw, nquad] = triangleIntegrationRule()
    
    for triIdx in range(meshdata.getTriangleCount()):
        detJ = meshdata.calcTriangleJacobianDeterminant(triIdx)
        Jinv = meshdata.calcTriangleInverseJacobian(triIdx)
        
        element_A = np.zeros((3*order,3*order))
        element_b = np.zeros((3*order,1))
        
        for i in range(element_A.shape[0]):
            for j in range (element_A.shape[1]):
                for k in range(nquad):
                    refpos = quadp[k][:]
                    temp = quadw[k]*np.abs(detJ)*np.dot(
                        # transpose this one to align shapes for dot product
                            np.dot(Jinv, delphi[i](*refpos)).transpose(),
                            np.dot(Jinv, delphi[j](*refpos))
                        ) 
                    element_A[i][j] = element_A[i][j] + temp
                    if j == 0:
                        mP = meshdata.calcTriangularIntegrationPoint(triIdx,quadp[k][:])
                        element_b[i] = element_b[i] + (
                            quadw[k]*forcingfunc(*mP)*np.abs(detJ)*phi[i](*refpos)
                            )
        nodes = meshdata.getTriangleNodes(triIdx,order)
        for i in range(element_A.shape[0]):
            g_i = nodes[i]
            for j in range(element_A.shape[0]):
                g_j = nodes[j]
                global_A[g_i][g_j] = global_A[g_i][g_j] + element_A[i][j]
            global_b[g_i] = global_b[g_i] + element_b[i]
    
    return global_A, global_b


def setBoundaryValues(meshdata, A, b, order):
    nodesD = set()
    nodesN = set()
    
    for lineIdx in range(meshdata.getLineCount()):
        tag = meshdata.getLineTag(lineIdx)
        if tag == 2:
            nodesD = nodesD.union(set(meshdata.getLineNodes(lineIdx, order)))
        elif tag == 3:
            nodesN = nodesN.union(set(meshdata.getLineNodes(lineIdx, order)))
        # else internal node or strange exotic boundary condition
    # size of b should be the number of nodes
    for i in range(b.shape[0]):
        if i in nodesD:
            # enforce dirichlet boundary
            # for this problem we can use lifting to make sure that they're zero
            b[i][0] = 0
            A[i][:] = 0
            A[i][i] = 1
        
        #do nothing for v.N. boundary conditions
    
    return A, b
            

mesh = MeshOperations("mesh/mesh1024.msh")
order = 1

if order > 2 or order < 1:
    order = 1

basis = [
       lambda xi, eta: 1-xi-eta,
       lambda xi, eta: xi,
       lambda xi, eta: eta
       ]
gradbasis = [
    lambda xi,eta: -1*np.ones((2,1)),
    lambda xi,eta: np.array([[1],[0]]),
    lambda xi,eta: np.array([[0],[1]])
    ]

def f(x,y):
    return 1

A, b = setEquationsSystem(mesh, basis, gradbasis, order, f)

A2, b2 = setBoundaryValues(mesh, A.copy(), b.copy(), order)

u = np.linalg.solve(A2, b2)

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot_trisurf(mesh.getListOfNodePositions()[:,0].flatten(),mesh.getListOfNodePositions()[:,1].flatten(),u.flatten())
plt.show()


    

