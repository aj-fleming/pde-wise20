# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 20:03:06 2020

@author: AJ
"""

import gmsh

gmsh.initialize()
gmsh.model.add("prog2")
geom = gmsh.model.geo

char_length = 0.1

# 2x2 square corners
geom.addPoint(0,0,0, char_length, 1)
geom.addPoint(2,0,0, char_length, 2)
geom.addPoint(2,2,0, char_length, 3)
geom.addPoint(0,2,0, char_length, 4)

# circle centered at (1,1)
geom.addPoint(1,1,0,char_length,5)
geom.addPoint(1,0.5,0,char_length,6)
geom.addPoint(1,1.5,0,char_length,7)

geom.addLine(1, 2, 1)
geom.addLine(2, 3, 2)
geom.addLine(3, 4, 3)
geom.addLine(4, 1, 4)

geom.addCircleArc(6, 5, 7, 5)
geom.addCircleArc(7, 5, 6, 6)

geom.addCurveLoop([1,2,3,4],1)
geom.addCurveLoop([5,6],2)

geom.addPlaneSurface([1,2],1)

geom.synchronize()

gmsh.model.addPhysicalGroup(1, [4], 2)
gmsh.model.addPhysicalGroup(1, [5,6], 4)
gmsh.model.addPhysicalGroup(1, [2], 5)
gmsh.model.addPhysicalGroup(1, [1,3], 3)
gmsh.model.addPhysicalGroup(2,[1], 1)

gmsh.fltk.run()

gmsh.finalize()