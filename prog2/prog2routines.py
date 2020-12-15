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

def setEquationsSystemQuadElements(meshdata, phi, delphi, forcingfunc):
    num_nodes = meshdata.getNodeCount()
    global_A = np.zeros((num_nodes, num_nodes))
    global_b = np.zeros((num_nodes, 1))
    
    qp1, qw, nq = triangleIntegrationRule()
    qp2 = np.add(1,-1*qp1)
    
    quads = meshdata.getListOfQuads()
    
    for quadIdx in range(meshdata.getQuadCount()):
        # define element A and b for first order quad elements 
        element_A = np.zeros((4,4))
        element_b = np.zeros((4,1))
        
        J1, J2 = meshdata.calcQuadJacobians(quadIdx)
        detJ1 = np.linalg.det(J1)
        detJ2 = np.linalg.det(J2)
        
        #find the inverses of the jacobian transposes
        J1_inv = np.linalg.inv(J1.transpose())
        J2_inv = np.linalg.inv(J2.transpose())
        nodes = quads[quadIdx]
        v = meshdata.getListOfNodePositions()[nodes]
        a1 = np.vstack([v[0,:]]*nq).transpose()
        a2 = np.vstack([v[1,:]-v[2,:]+v[3,:]]*nq).transpose()
        
        mp1 = (a1 + np.dot(J1,qp1.transpose())).transpose()
        mp2 = (a2 + np.dot(J2,qp2.transpose())).transpose()
        for i in range(4):
            for j in range(4):
                
                for k in range(nq):
                    refpos = qp1[k][:]
                    temp1 = qw[k]*np.abs(detJ1)*np.dot(
                        # transpose this one to align shapes for dot product
                            np.dot(J1_inv, delphi[i](*refpos)).transpose(),
                            np.dot(J1_inv, delphi[j](*refpos))
                        )
                    refpos2 = qp2[k,:]
                    temp2 = qw[k]*np.abs(detJ2)*np.dot(
                        # transpose this one to align shapes for dot product
                            np.dot(J2_inv, delphi[i](*refpos2)).transpose(),
                            np.dot(J2_inv, delphi[j](*refpos2))
                        )
                    element_A[i][j] = element_A[i][j] + temp1 + temp2
                    if j == 0:
                        temp1 = qw[k]*forcingfunc(*mp1[k])*np.abs(detJ1)*phi[i](*refpos)
                        temp2 = qw[k]*forcingfunc(*mp2[k])*np.abs(detJ2)*phi[i](*refpos2)
                        element_b[i] = element_b[i] + temp1 + temp2
        for i in range(4):
            g_i = nodes[i]
            for j in range(4):
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