# -*- coding: utf-8 -*-

from meshops import MeshOperations, triangleIntegrationRule
import numpy as np

def setEquationsSystem(meshdata, phi, delphi, order):
    global_A = np.zeros(meshdata.getNodeCount())
    global_b = np.zeros((meshdata.getNodeCount(),1))
    
    [quadp, quadw, nquad] = triangleIntegrationRule()
    
    for triIdx in range(meshdata.getTriangleCount()):
        detJ = meshdata.calcTriangleJacobianDeterminant(triIdx)
        Jinv = meshdata.calcTriangleInverseJacobian(triIdx)
        
        element_A = np.zeros(3*order)
        element_b = np.zeros(3*order)
        
        for i in range(element_A.shape[0]):
            for j in range (element_A.shape[1]):
                for k in range(nquad):
                    ref_x = quadp[k][0]
                    ref_y = quadp[k][1]
                    element_A[i][j] = element_A[i][j] + quadw[k]*(
                        np.dot(Jinv, delphi[i](ref_x, ref_y))
                        )
                    
def setBoundaryValues(mops):
    pass

mesh = MeshOperations("mesh/unitSquare1.msh")
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

setEquationsSystem(mesh, basis, gradbasis, order)



    

