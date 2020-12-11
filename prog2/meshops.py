# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 14:36:02 2020

@author: AJ Fleming
"""
import numpy as np

NODES_PER_ELEMENT_TYPE = [2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13]

def lineIntegrationRule():
    qp = np.zeros((3,1))
    qp[0] = -np.sqrt(3/5);
    qp[2] = np.sqrt(3/5);
    
    qw = np.zeros((3,1))
    qw[0] = 5/9
    qw[1] = 8/9
    qw[2] = 5/9
    
    return qp, qw, 3

def triangleIntegrationRule():
    """
    Returns
    -------
    quadPoints : 7x2 ndarray of floats
        quadrature points for the reference triangle stored as row vectors
    quadWeights : 7x1 array of floats
        weights for the given quadrature points
    int
        7
    """
    quadPoints = np.zeros((7,2))
    # points in the xi direction
    quadPoints[0][0] = 1/3;
    quadPoints[1][0] = (6 + np.sqrt(15))/21;
    quadPoints[2][0] = (9 - 2*np.sqrt(15))/21;
    quadPoints[3][0] = (6 + np.sqrt(15))/21;
    quadPoints[4][0] = (6 - np.sqrt(15))/21;
    quadPoints[5][0] = (9 + 2*np.sqrt(15))/21;
    quadPoints[6][0] = (6 - np.sqrt(15))/21;
    # eta-values
    quadPoints[0][1] = 1/3;
    quadPoints[1][1] = (6 + np.sqrt(15))/21;
    quadPoints[2][1] = (6 + np.sqrt(15))/21;
    quadPoints[3][1] = (9 - 2*np.sqrt(15))/21;
    quadPoints[4][1] = (6 - np.sqrt(15))/21;
    quadPoints[5][1] = (6 - np.sqrt(15))/21;
    quadPoints[6][1] = (9 + 2*np.sqrt(15))/21;
    # weights
    quadWeights = np.zeros((7,1));
    quadWeights[0] = 9/80;
    quadWeights[1:4] = (155 + np.sqrt(15))/2400;
    quadWeights[4:7] = (155 - np.sqrt(15))/2400;
    
    return quadPoints, quadWeights, 7


class MeshOperations:
    def __init__(self, meshfile):
        """
        Loads the mesh from disk. Contains many useful methods for working with
        2-d triangular meshes.

        Parameters
        ----------
        meshfile : file path
            location of an ASCII gmesh2 file.

        Returns
        -------
        None.

        """
        self.mesh = MeshGrid(meshfile)          
    
    def reset(self):
        self.mesh.reload()
    
    def __repr__(self):
        bbinf = "Mesh bounded by ({:4f},{:4f},{:4f}) and ({:4f},{:4f},{:4f}), with:".format(*self.mesh.bounding_box[0],
                                                                                            *self.mesh.bounding_box[1])
        nodeinf = "\t{} total nodes".format(self.mesh.num_nodes)
        lineinf ="\t{} total lines and {} second order lines".format(self.mesh.num_lines,self.mesh.num_lines3)
        triinf = "\t{} total triangles and {} second order triangles".format(self.mesh.num_tris, self.mesh.num_tris6)
        return "\n".join([bbinf, nodeinf, lineinf, triinf])
    
    def getLineTag(self, i):
        """
        Parameters
        ----------
        i : integer
            the edge element whose tag we want. indexing starts at 0!

        Returns
        -------
        integer
            1 for inner edges
            2 for outer edge with a dirichlet boundary value
            3 for an outer edge with a von Neumann boundary value

        """
        return self.mesh.lines[i][2]
    
    def getTriangleTag(self, i):
        return self.mesh.tris[i][3]
    
    def getNodeCount(self):
        return self.mesh.num_nodes
    
    def getListOfTriangles(self):
        return self.mesh.tris[0:self.mesh.num_tris][0:2]
    
    def getListOfNodePositions(self):
        return self.mesh.node_positions[:][0:2]
    
    def getTriangleCount(self):
        return self.mesh.num_tris
    
    def getLineCount(self):
        return self.mesh.num_lines
    
    def getLineNodes(self, lineIdx, order):
        if order == 1:
            return self.mesh.lines[lineIdx][0:2]
        return self.mesh.lines3[lineIdx][0:3]
    
    def getTriangleNodes(self, triIdx, order):
        if order == 1:
            return self.mesh.tris[triIdx][0:3]
        return self.mesh.tris6[triIdx][0:4]
    
    def calcLineJacobian(self, lineIdx):
        endpoints = self.mesh.lines[lineIdx][0:2]
        points = self.mesh.node_positions[endpoints][0:2]
        
        J = np.zeros((1,2))
        J[0] = 0.5 * (points[1][0] - points[0][0])
        J[1] = 0.5 * (points[1][1] - points[0][1])
        return J
    
    def calcTriangleJacobian(self, triIdx):
        vertexIds = self.mesh.tris[triIdx][0:3]
        points = self.mesh.node_positions[vertexIds][0:3] 
        
        J = np.zeros((2,2))
        J[0][0] = points[1][0] - points[0][0]
        J[0][1] = points[1][1] - points[0][1]
        J[1][0] = points[2][0] - points[0][0]
        J[1][1] = points[2][1] - points[0][1]
        
        return J
    
    def calcTriangleInverseJacobian(self, triIdx):
        return np.linalg.inv(self.calcTriangleJacobian(triIdx))
    
    def calcLineJacobianDeterminant(self, lineIdx):
        J = self.calcLineJacobian(lineIdx)
        return np.sqrt(np.sum(np.power(J,2)))
    
    def calcTriangleJacobianDeterminant(self, triIdx):
        return np.abs(np.linalg.det(self.calcTriangleJacobian(triIdx)))
    
    def calcTriangularIntegrationPoints(self, triIdx, int_point):
        J = self.calcTriangleJacobian(triIdx)
        p1n = self.mesh.tris[triIdx][0]
        refPos = self.mesh.node_positions[p1n][0:2]
        
        tv = np.dot(J.transpose(), int_point)
        return np.add(refPos, tv)
    
    
class MeshGrid:
    """
    Describes a GMesh mesh. This is a low-level structure to keep track of small details about the mesh.
    Maintins element tag and connectivity information.
    """
    def __init__(self, file):
        self.filename = file
        self.bounding_box = np.zeros((2,3))
        self.initialized = self.reload()

    def reload(self):
        """
        Load the mesh file from disk (perhaps again if it has been modified in the program)

        Returns
        -------
        bool
            True if the v2 ASCII mesh file was successfully loaded.

        """
        with open(self.filename,'r') as meshfile:
            # scan file until we reach a mesh format declarator
            if not scan_for_keyword(meshfile, "$meshformat"):
                return False
            print("found format")
            # read mesh format information
            self.meshformat = meshfile.readline()
            #check for end of mesh formatting block
            if meshfile.readline().lower().strip() != "$endmeshformat":
                print("Can only read ASCII meshes.")
                return False

            if not scan_for_keyword(meshfile, "$nodes"):
                return False
            print("reading nodes")

            self.num_nodes = int(meshfile.readline())
            self.node_positions = np.zeros((self.num_nodes, 3))
            nodeids = [0]*self.num_nodes
            for i in range(self.num_nodes):
                nodeinf = meshfile.readline().split()
                # shift to zero-indexing from gmsh/matlab 1-indexing
                nodeids[i] = int(nodeinf[0]) - 1
                nodex = np.array([float(k) for k in nodeinf[1:]])
                #set axis-aligned bounding box for the mesh
                if (i == 0):
                    self.bounding_box[0] = nodex
                    self.bounding_box[1] = nodex
                else:
                    self.bounding_box[0] = [min(self.bounding_box[0][k],nodex[k]) for k in range(3)]
                    self.bounding_box[1] = [max(self.bounding_box[1][k],nodex[k]) for k in range(3)]
                self.node_positions[i] = nodex
            print("read {} nodes".format(self.num_nodes))
            if not scan_for_keyword(meshfile, "$endnodes"):
                return False
            if not scan_for_keyword(meshfile, "$elements"):
                return False
            print("found elements")

            self.num_elements = int(meshfile.readline())
            #constants given by the file format
            num_infos = 4
            tagidx = 3
            self.element_infos = [[0]*num_infos]*self.num_elements
            self.element_tags = [0]*self.num_elements
            self.num_points = 0
            self.num_lines = 0
            self.num_tris = 0
            self.num_quads = 0
            # self.num_tets = 0
            # self.num_hexas = 0
            # self.num_prisms = 0
            # self.num_pyramids = 0
            self.num_lines3 = 0
            self.num_tris6 = 0

            self.points = np.zeros((self.num_elements,2), np.int32)
            self.lines = np.zeros((self.num_elements,3), np.int32)
            self.tris = np.zeros((self.num_elements,4), np.int32)
            self.quads = np.zeros((self.num_elements,5), np.int32)
            # self.tets = np.zeros((self.num_elements,5), np.int32)
            # self.hexas = np.zeros((self.num_elements,9), np.int32)
            # self.prisms = np.zeros((self.num_elements,7), np.int32)
            # self.pyramids = np.zeros((self.num_elements,6), np.int32)
            self.lines3 = np.zeros((self.num_elements,4), np.int32)
            self.tris6 = np.zeros((self.num_elements,7), np.int32)

            tokens = []
            tline = meshfile.readline().lower().strip()
            while tline != "$endelements":
                if not tline:
                    return False
                tokens = tokens + [int(k) for k in tline.split()]
                tline = meshfile.readline().lower().strip()
            for i in range(self.num_elements):
                self.element_infos[i] = [tokens.pop(0) for k in range(num_infos)]
                # I have honestly no clue what this means, but it consumes tokens
                #   so it's staying in the code
                self.element_tags[i] = [tokens.pop(0) for k in range(self.element_infos[i][2]-1)]
                # minus 1s to shift from one-indexing to zero-indexing
                element_nodes = [tokens.pop(0)-1 for k in range(NODES_PER_ELEMENT_TYPE[self.element_infos[i][1]-1])]

                if self.element_infos[i][1] == 15:
                    self.points[self.num_points][0] = nodeids[element_nodes[0]]
                    self.points[self.num_points][1] = self.element_infos[i][tagidx]
                    self.num_points = self.num_points + 1
                elif self.element_infos[i][1] == 1:
                    self.add_line(i, nodeids, element_nodes, 1)
                elif self.element_infos[i][1] == 8:
                    self.add_line(i, nodeids, element_nodes, 2)
                elif self.element_infos[i][1] == 2:
                    self.add_triangle(i, nodeids, element_nodes, 1)
                elif self.element_infos[i][1] == 9:
                    self.add_triangle(i, nodeids, element_nodes, 2)
                elif self.element_infos[i][1] == 3:
                    for j in range(4):
                        self.quads[self.num_quads][j] = nodeids[element_nodes[j]]
                    self.quads[self.num_quads][3] = self.element_infos[i][tagidx]
                    self.num_quads = self.num_quads + 1

                #TODO tetras/hexes/prisms/pyramids
                

        return True
        
    def add_line(self, i, nodeids, element_nodes, order):
        self.lines[self.num_lines][0] = nodeids[element_nodes[0]]
        self.lines[self.num_lines][1] = nodeids[element_nodes[1]]
        self.lines[self.num_lines][2] = self.element_infos[i][3]
        self.num_lines = self.num_lines + 1
        if order == 2:
            self.lines3[self.num_lines3][0] = nodeids[element_nodes[0]]
            self.lines3[self.num_lines3][1] = nodeids[element_nodes[1]]
            self.lines3[self.num_lines3][2] = nodeids[element_nodes[2]]
            self.lines3[self.num_lines3][3] = self.element_infos[i][3]
            self.num_lines3 = self.num_lines3 + 1
   
    
    def add_triangle(self, i, nodeids, element_nodes, order):
        for j in range(3):
            self.tris[self.num_tris][j] = nodeids[element_nodes[j]]
        self.tris[self.num_tris][3] = self.element_infos[i][3]
        self.num_tris = self.num_tris + 1
        if order == 2:
            for j in range(6):
                self.tris6[self.num_tris6][j] = nodeids[element_nodes[j]]
            self.tris6[self.num_tris6][6] = self.element_infos[i][3]
            self.num_tris6 = self.num_tris6 + 1    

def scan_for_keyword(file, keyword):
    tline = file.readline().lower().strip()
    while tline != keyword.lower().strip():
        if not tline:
            return False
        tline = file.readline().lower().strip()
    return True
