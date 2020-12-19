# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 20:40:19 2020

@author: AJ
"""

import gmsh

gmsh.initialize()
gmsh.model.add("prog2")

char_length = 1

gmsh.model.geo.addPoint(0,0,0, char_length, 1)
gmsh.model.geo.addPoint(1,0,0, char_length, 2)
gmsh.model.geo.addPoint(1,1,0, char_length, 3)
gmsh.model.geo.addPoint(0,1,0, char_length, 4)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addCurveLoop([1,2,3,4], 1)

gmsh.model.geo.addPlaneSurface([1],1)




# for i in range(4):
#     gmsh.model.geo.mesh.setTransfiniteCurve(i+1, 10)
# gmsh.model.geo.mesh.setTransfiniteSurface(1,cornerTags=[1,2,3,4])
gmsh.model.geo.mesh.setRecombine(2,1)

i_tag = 1
d_tag = 2
vn_tag = 3

gmsh.model.geo.synchronize()


gmsh.model.addPhysicalGroup(1,[1,3,4],d_tag)
gmsh.model.addPhysicalGroup(1,[2],vn_tag)
gmsh.model.addPhysicalGroup(2,[1], 1)


gmsh.fltk.run()

gmsh.finalize()

