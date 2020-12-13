# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 14:58:02 2020

@author: AJ
"""
import numpy as np
from meshops import triangleIntegrationRule


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