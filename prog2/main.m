%% Preparation: add some paths
mops = MeshOperations('mesh/unitSquare1.msh');
mesh_A = MeshOperations('mesh/unitSquare1.msh');
mesh_B = MeshOperations('mesh/unitSquare2.msh');
order = 1;

if order > 2 || order < 1
    order = 1;
end

basis1stOrder = [
       lambda xi, eta: 1-xi-eta,
       lambda xi, eta: xi,
       lambda xi, eta: eta
       ]
basis2ndOrder = [];
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

A1a2, b1a2 = setBoundaryValues(mesh_A, A1a.copy(), b1a.copy(), order);
A1b2, b1b2 = setBoundaryValues(mesh_B, A1b.copy(), b1b.copy(), order);
A1c2, b1c2 = setBoundaryValues(mesh_B, A1c.copy(), b1c.copy(), order);

u1a = np.linalg.solve(A1a2, b1a2)
u1b = np.linalg.solve(A1b2, b1b2)
u1c = np.linalg.solve(A1c2, b1c2)

fig = plt.figure(1)
fig.suptitle("Solution to 1a: f(x,y) = 1 on unitSquare1")
mesh_A.plot(fig.gca(projection='3d'), u1a)
fig = plt.figure(2)
fig.suptitle("Solution to 1b: f(x,y) = 1 on unitSquare2")
mesh_B.plot(fig.gca(projection='3d'), u1b)
fig = plt.figure(3)
fig.suptitle("Solution to 1c: f(x,y) = sin(2*pi*x) on unitSquare2")
mesh_B.plot(fig.gca(projection='3d'), u1c)

mesh64 = MeshOperations("mesh/mesh64.msh")
mesh256 = MeshOperations("mesh/mesh256.msh")
mesh1024 = MeshOperations("mesh/mesh1024.msh")

a64, b64 = setBoundaryValues(mesh64,
                             *setEquationsSystem(mesh64,basis1stOrder,gradBasis1stOrder,order,f1a),
                             order)
a256, b256 = setBoundaryValues(mesh256,
                             *setEquationsSystem(mesh256,basis1stOrder,gradBasis1stOrder,order,f1a),
                             order)
a1024, b1024 = setBoundaryValues(mesh1024,
                             *setEquationsSystem(mesh1024,basis1stOrder,gradBasis1stOrder,order,f1a),
                             order)
u64=a64\b64;
u256 = a256\b256;
u1024 = a1024\b1024;

fig = plt.figure(4)
fig.suptitle("Solution to $\Delta u = 1$ on mesh64")
mesh64.plot(fig.gca(projection='3d'), u64)
fig = plt.figure(5)
fig.suptitle("Solution to $\Delta u = 1$ on mesh256")
mesh256.plot(fig.gca(projection='3d'), u256)
fig = plt.figure(6)
fig.suptitle("Solution to $\Delta u = 1$ on mesh1024")
mesh1024.plot(fig.gca(projection='3d'), u1024)

errs = [
        mesh64.calcL2ErrorPoisson(u64,order,basis1stOrder),
        mesh256.calcL2ErrorPoisson(u256,1,basis1stOrder),
        mesh1024.calcL2ErrorPoisson(u1024,1,basis1stOrder)
        ]
plt.show()
function [A, b] = setEquationsSystem(meshdata, phi, delphi, order, forcingfunc)
    A = np.zeros(meshdata.getNodeCount(),meshdata.getNodeCount());
    b = np.zeros(meshdata.getNodeCount(),1);
    
    [quadp, quadw, nquad] = triangleIntegrationRule()
    
    for triIdx=1:meshdata.getTriangleCount()
        detJ = meshdata.calcTriangleJacobianDeterminant(triIdx)
        Jinv = meshdata.calcTriangleInverseJacobian(triIdx)
        
        element_A = np.zeros(3*order,3*order);
        element_b = np.zeros(3*order,1);
        
        for i=1:element_A.shape(1)
            for j=1:element_A.shape(2)
                for k=1:nquad
                    refpos = quadp(k,:);
                    temp = quadw(k)*abs(detJ)*dot(transpose(dot(Jinv, delphi(i)*refpos)),dot(Jinv, delphi(j)*refpos));
                    element_A(i,j) = element_A(i,j) + temp;
                    if j == 0
                        mP = meshdata.calcTriangularIntegrationPoint(triIdx,quadp(k,:));
                        element_b(i) = element_b(i) + (quadw(k)*forcingfunc*mP*abs(detJ)*phi(i)*refpos);
                    end
                end
            end
        end
        nodes = meshdata.getTriangleNodes(triIdx,order);
        for i=1:element_A.shape(0)
            g_i = nodes(i);
            for j=1:element_A.shape(0)
                g_j = nodes(j);
                global_A(g_i,g_j) = global_A(g_i,g_j) + element_A(i,j);
            global_b(g_i) = global_b(g_i) + element_b(i);
            end
        end
    end
end


function [A,b] = setBoundaryValues(meshdata, A, b, order)
    nodesD = [];
    nodesN = [];
    
    for lineIdx=1:meshdata.getLineCount()
        tag = meshdata.getLineTag(lineIdx);
        if tag == 2
            nodesD = union(nodesD,meshdata.getLineNodes(lineIdx, order));
        elseif tag == 3
            nodesN = union(nodesN,meshdata.getLineNodes(lineIdx, order));
        end
    end
%         # else internal node or strange exotic boundary condition
%     # size of b should be the number of nodes
    for i=1:b.shape(0)
        if ismember(i,nodesD)
%             # enforce dirichlet boundary
%             # for this problem we can use lifting to make sure that they're zero
            b(i,0) = 0;
            A(i,:) = 0;
            A(i,i) = 1;
            
        end
    end
 
    
end