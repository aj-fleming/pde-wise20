# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 21:39:00 2020

@author: AJ
"""

import matplotlib.pyplot as plt
import numpy as np

from meshops import triangleIntegrationRule

# get the integration points for the triangle on the lower left
qP, qW, n = triangleIntegrationRule()
# create integration points for the triangle on the upper right
tqP = np.add(1,-1*qP)
fig = plt.figure()
ax = fig.gca()
ax.set_aspect('equal')

ax.scatter(qP[:,0],qP[:,1], marker="o", color="red")
ax.scatter(tqP[:,0],tqP[:,1], marker="o", color="blue")
ax.plot([0,0,1],[1,0,0], color="red")
ax.plot([1,1,0],[0,1,1],color="blue")
ax.plot([1,0],[0,1],color="purple",linestyle=":")

# vertices of arbitrary quad
v = np.array([
    [2,2],
    [3,1],
    [4,3],
    [2,4]
    ])

#vertices of unit square
v2 = np.array([
    [0,0],
    [1,0],
    [1,1],
    [0,1]
    ])
J1 = np.zeros((2,2))
J2 = np.zeros((2,2))
        
J1[0][0] = v[1][0] - v[0][0]
J1[1][0] = v[1][1] - v[0][1]
J1[0][1] = v[3][0] - v[0][0]
J1[1][1] = v[3][1] - v[0][1]

J2[0][0] = v[2][0] - v[3][0]
J2[0][1] = v[2][0] - v[1][0]
J2[1][0] = v[2][1] - v[3][1]
J2[1][1] = v[2][1] - v[1][1]
        
a = v[1,:]-v[2,:]+v[3,:]
A = np.vstack([a]*n).transpose()
points = A + np.dot(J2,tqP.transpose())
verts = [1,2,3]
A = np.vstack([a]*3).transpose()
boundary = A + np.dot(J2, v2[verts].transpose())
ax.scatter(points[0,:],points[1,:],color="blue", marker="+")
ax.plot(boundary[0,:],boundary[1,:],color="blue")

a = v[0,:]
A = np.vstack([a]*n).transpose()
points = A + np.dot(J1,qP.transpose())
verts = [3,0,1]
A = np.vstack([a]*3).transpose()
boundary = A + np.dot(J1, v2[verts].transpose())
ax.scatter(points[0,:],points[1,:],color="red",marker="+")
ax.plot(boundary[0,:],boundary[1,:],color="red")
ax.plot(boundary[0,[0,2]], boundary[1,[0,2]],color="purple", linestyle=":")



plt.show()