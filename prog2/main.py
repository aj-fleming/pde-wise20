# -*- coding: utf-8 -*-

from meshops import MeshOperations, triangleIntegrationRule
import numpy as np
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
            

mesh_A = MeshOperations("mesh/unitSquare1.msh")
mesh_B = MeshOperations("mesh/unitSquare2.msh")
order = 1

if order > 2 or order < 1:
    order = 1

basis1stOrder = [
       lambda xi, eta: 1-xi-eta,
       lambda xi, eta: xi,
       lambda xi, eta: eta
       ]
basis2ndOrder = [
       
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

A1a, b1a = setEquationsSystem(mesh_A, basis1stOrder, gradBasis1stOrder, order, f1a)
A1b, b1b = setEquationsSystem(mesh_B, basis1stOrder, gradBasis1stOrder, order, f1a)
A1c, b1c = setEquationsSystem(mesh_B, basis1stOrder, gradBasis1stOrder, order, f1c)

A1a2, b1a2 = setBoundaryValues(mesh_A, A1a.copy(), b1a.copy(), order)
A1b2, b1b2 = setBoundaryValues(mesh_B, A1b.copy(), b1b.copy(), order)
A1c2, b1c2 = setBoundaryValues(mesh_B, A1c.copy(), b1c.copy(), order)

u1a = np.linalg.solve(A1a2, b1a2)
u1b = np.linalg.solve(A1b2, b1b2)
u1c = np.linalg.solve(A1c2, b1c2)

fig1 = plt.figure(1)
mesh_A.plot(fig1.gca(projection='3d'), u1a)
fig2 = plt.figure(2)
mesh_B.plot(fig2.gca(projection='3d'), u1b)
fig3 = plt.figure(3)
mesh_B.plot(fig3.gca(projection='3d'), u1c)
plt.show()

print(mesh_A.calcL2ErrorPoisson(u1a, 1, basis1stOrder))


    

