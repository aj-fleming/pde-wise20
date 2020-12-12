classdef MeshOperations < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %MeshOperations Mesh Wrapper Class
    %   Contains the struct with mesh information by getting the mesh path
    %   and provides element transformation and integration rules
    %   for specific element geometries
    %   All procedures that are specific to the element geometry but not
    %   to the shape functions are gathered here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties(GetAccess = private)
        %has Mesh Data from gmsh interface
        mesh;
    end
    
    methods
        function meshOperationsObj = MeshOperations(meshPath)
            meshOperationsObj.mesh = meshOperationsObj.load_gmsh(meshPath);
            %creating low order data out of second order mesh
            if(meshOperationsObj.mesh.nbTriangles6 > 0)
                meshOperationsObj.mesh.nbTriangles = meshOperationsObj.mesh.nbTriangles6;
                meshOperationsObj.mesh.nbLines = meshOperationsObj.mesh.nbLines3;
                %
                meshOperationsObj.mesh.TRIANGLES = zeros(meshOperationsObj.mesh.nbTriangles6, 7);
                %first three points of triangle6 are low order points
                %(corners of the triangle)
                meshOperationsObj.mesh.TRIANGLES(1:meshOperationsObj.mesh.nbTriangles6, 1:3) = meshOperationsObj.mesh.TRIANGLES6(1:meshOperationsObj.mesh.nbTriangles6, 1:3);
                meshOperationsObj.mesh.TRIANGLES(1:meshOperationsObj.mesh.nbTriangles6, 4) = meshOperationsObj.mesh.TRIANGLES6(1:meshOperationsObj.mesh.nbTriangles6, 7);
                %
                meshOperationsObj.mesh.LINES = zeros(meshOperationsObj.mesh.nbLines3, 4);
                %first two points of triangle3 are low order points
                %(bounds of the line)
                meshOperationsObj.mesh.LINES(1:meshOperationsObj.mesh.nbLines3, 1:2) = meshOperationsObj.mesh.LINES3(1:meshOperationsObj.mesh.nbLines3, 1:2);
                meshOperationsObj.mesh.LINES(1:meshOperationsObj.mesh.nbLines3, 3) = meshOperationsObj.mesh.LINES3(1:meshOperationsObj.mesh.nbLines3, 4);
            end
        end
        
        function getMeshInformation(meshOperationsObj)
            fprintf('################################################ \n')
            fprintf('## number of volume elements (triangles): %d \n', meshOperationsObj.mesh.nbTriangles)
            fprintf('## number of boundary elements (lines): %d \n', meshOperationsObj.mesh.nbLines)
            fprintf('################################################ \n')
        end
        
        function reset(meshOperationsObj, meshPath)
            meshOperationsObj.mesh = meshOperationsObj.load_gmsh(meshPath);
        end
        
        function plot(meshOperationsObj, u, component)
            triangles = meshOperationsObj.getVolumeElementListTriangles;
            nodes = meshOperationsObj.getNodeList();
            nbNodes = meshOperationsObj.getNumberNodes();
            trisurf(triangles(:,1:3), nodes(:,1), nodes(:,2), u((component-1)*nbNodes+1:component*nbNodes));
             title( ['Solution u_' num2str( component )] )
            xlabel('x');
            ylabel('y');
            set(gca,'FontSize',20);
        end
        
        function numElems = getNumberOfLines(meshOperationsObj)
            numElems = meshOperationsObj.mesh.nbLines;
        end
        
        function numElems = getNumberOfTriangles(meshOperationsObj)
            numElems = meshOperationsObj.mesh.nbTriangles;
        end
        
        function tag = getTagOfLine(meshOperationsObj, elementNumber)
            %delivers "physical" tag
            tag = meshOperationsObj.mesh.LINES(elementNumber, 3);
        end
        function tag = getTagOfTriangle(meshOperationsObj, elementNumber)
            %delivers "physical" tag
            tag = meshOperationsObj.mesh.TRIANGLES(elementNumber, 4);
        end
        
        function numNodes = getNumberNodes(meshOperationsObj)
            numNodes = meshOperationsObj.mesh.nbNod;
        end
        
        function volElems = getVolumeElementListTriangles(meshOperationsObj)
            %number of elements
            nbElems = meshOperationsObj.getNumberOfTriangles;
            %tags are not included
            volElems = meshOperationsObj.mesh.TRIANGLES(1:nbElems, 1:3);
        end
        
        function nodes = getNodeList(meshOperationsObj)
            %only 2D-Nodes
            nodes = meshOperationsObj.mesh.POS(:,1:2);
        end
        
        function nodeData = getNodeNumbersOfLine(meshOperationsObj, elementNumber, order)
            if(order == 1)
                nodeData = meshOperationsObj.mesh.LINES(elementNumber,1:2);
            elseif(order == 2)
                nodeData = meshOperationsObj.mesh.LINES3(elementNumber,1:3);
            end
        end
        
        function nodeData = getNodeNumbersOfTriangle(meshOperationsObj, elementNumber, order)
            if(order == 1)
                nodeData = meshOperationsObj.mesh.TRIANGLES(elementNumber,1:3);
            elseif(order == 2)
                nodeData = meshOperationsObj.mesh.TRIANGLES6(elementNumber,1:6);
            end
        end
        
        function jacobian = calcJacobianOfLine(meshOperationsObj, elementNumber)
            %identify nodes belonging to the line segment
            pointNumbers = meshOperationsObj.mesh.LINES(elementNumber,1:2);
            points = zeros(2,2);
            %reading x- and y-coordinates
            points(:,:) = meshOperationsObj.mesh.POS(pointNumbers,1:2);
            %Compute entries of Jacobian (just Vector!!!)
            jacobian = zeros(2,1);
            jacobian(1) = 0.5*(points(2,1) - points(1,1)); %dxdxi 1/2 * (x2-x1)
            jacobian(2) = 0.5*(points(2,2) - points(1,2)); %dydxi 1/2 * (y2-y1)
        end
        
        function jacobian = calcJacobianOfTriangle(meshOperationsObj, elementNumber)
            %identify nodes belonging to the triangle
            pointNumbers = meshOperationsObj.mesh.TRIANGLES(elementNumber,1:3);
            points = zeros(3,2);
            %reading x- and y-coordinates
            points(:,:) = meshOperationsObj.mesh.POS(pointNumbers,1:2);
            %Compute entries of Jacobian
            jacobian = zeros(2,2);
            jacobian(1,1) = points(2,1) - points(1,1); %dxdxi (x2-x1)
            jacobian(1,2) = points(2,2) - points(1,2); %dydxi (y2-y1)
            jacobian(2,1) = points(3,1) - points(1,1); %dxdeta (x3-x1)
            jacobian(2,2) = points(3,2) - points(1,2); %dydeta (y3-y1)
        end
        
        function inverseJacobian = calcInverseJacobianOfTriangle(meshOperationsObj, elementNumber)
            jacMat = meshOperationsObj.calcJacobianOfTriangle(elementNumber);
            inverseJacobian = inv(jacMat);
        end
        
        function jacobianDet = calcJacobianDeterminantOfLine(meshOperationsObj, elementNumber)
            jacMat = meshOperationsObj.calcJacobianOfLine(elementNumber);
            jacobianDet = sqrt((jacMat(1,1))^2 + (jacMat(2,1))^2);
        end
        
        function jacobianDet = calcJacobianDeterminantOfTriangle(meshOperationsObj, elementNumber)
            jacMat = meshOperationsObj.calcJacobianOfTriangle(elementNumber);
            jacobianDet = abs(det(jacMat));
        end
        
        function mappedIP = calcMappedIntegrationPointOfLine(meshOperationsObj, elementNumber, IP)
            jacMat = meshOperationsObj.calcJacobianOfLine(elementNumber);
            pointNumbers = meshOperationsObj.mesh.LINES(elementNumber,1:2);
            points = zeros(2,2);
            %reading x- and y-coordinates
            points(:,:) = meshOperationsObj.mesh.POS(pointNumbers,1:2);
            %different reference position due to non unit reference
            %interval
            mappedIP = IP(1) * jacMat + 0.5*[points(2,1) + points(1,1); points(2,2) + points(1,2)];
        end
        
        function mappedIP = calcMappedIntegrationPointOfTriangle(meshOperationsObj, elementNumber, IP)
            jacMat = meshOperationsObj.calcJacobianOfTriangle(elementNumber);
            point1Number = meshOperationsObj.mesh.TRIANGLES(elementNumber,1);
            refPos = meshOperationsObj.mesh.POS(point1Number,1:2);
            transMat =  transpose(jacMat);
            tempVec(1) =  dot(transMat(1,:),IP);
            tempVec(2) =  dot(transMat(2,:),IP);
            mappedIP = tempVec + refPos;
        end
        
        %look out for orientation during geometry setup
        %counter clockwise for outward facing normal
        %clockwise for inside facing normal
        function normalVec = getNormalVectorOfLine(meshOperationsObj, elementNumber)
            pointNumbers = meshOperationsObj.mesh.LINES(elementNumber,1:2);
            points = zeros(2,2);
            %reading x- and y-coordinates
            points(:,:) = meshOperationsObj.mesh.POS(pointNumbers,1:2);
            %check direction here
            tangentVec = points(2,:) - points(1,:);
            magnitude = sqrt((tangentVec(1))^2 + (tangentVec(2))^2);
            normalVec(1) = tangentVec(2)/magnitude;
            normalVec(2) = -tangentVec(1)/magnitude;
        end
        
        function l2Error = L2errorOfPoissonProblem(meshOperationsObj, u, order)
            nTriangles = meshOperationsObj.getNumberOfTriangles();
            l2Error = 0;
            
            % Integration rule
            [quadWeights, quadPoints, numIPs] = meshOperationsObj.IntegrationRuleOfTriangle();
            
            % Loop over each triangle
            for i = 1:nTriangles
                % Loop over element integration points
                determinant_of_Jacobian = meshOperationsObj.calcJacobianDeterminantOfTriangle(i);
                for j = 1:numIPs
                    if (order == 1)
                        shape = zeros(3,1);
                        shape(1) = 1 - quadPoints(j,1) - quadPoints(j,2);
                        shape(2) = quadPoints(j,1);
                        shape(3) = quadPoints(j,2);
                    elseif (order == 2)
                        shape = zeros(6,1);
                        shape(1) = (1 - quadPoints(j,1) - quadPoints(j,2))*(1 - 2*quadPoints(j,1) - 2*quadPoints(j,2));
                        shape(2) = quadPoints(j,1)*(2*quadPoints(j,1) - 1);
                        shape(3) = quadPoints(j,2)*(2*quadPoints(j,2) - 1);
                        shape(4) = 4*quadPoints(j,1)*(1 - quadPoints(j,1) - quadPoints(j,2));
                        shape(5) = 4*quadPoints(j,1)*quadPoints(j,2);
                        shape(6) = 4*quadPoints(j,2)*(1 - quadPoints(j,1) - quadPoints(j,2));
                    end
                    % Calculate analytical value
                    mappedIP = meshOperationsObj.calcMappedIntegrationPointOfTriangle(i, quadPoints(j,:));
                    u_exact = 0;
                    fourierOrder = 30;
                    for k = 1:fourierOrder
                        for l = 1:fourierOrder
                            if(mod(l,2) ~= 0)                                                                                                       % coefficient and 
                                coeff = 16/((pi)^4) * 1 / (l^3*(2*k-1) + l*(2*k-1)^3/4);                                        % basisfunction for
                                u_exact = u_exact + coeff * sin((k-0.5)*pi*mappedIP(1)).*sin(l*pi*mappedIP(2));  % neumann bc (=0) at one face
%                             if((mod(k,2) ~= 0) && (mod(l,2) ~= 0))                                                                     % coefficient and
%                                 coeff = 16/((k^2+l^2)*k*l*(pi)^4);                                                                       % basisfunction for
%                                 u_exact = u_exact + coeff * sin(k*pi*mappedIP(1))*sin(l*pi*mappedIP(2));            % dirichlet bc (=0) on all faces
                            end
                        end
                    end
                    
                    % local L2 Error
                    conVector = meshOperationsObj.getNodeNumbersOfTriangle(i, order);
                    errorVal = (dot(shape, u(conVector)) - u_exact)^2;
                    l2Error = l2Error + quadWeights(j) * determinant_of_Jacobian * errorVal;
                end
            end
            l2Error = sqrt(l2Error);
        end
    end
    
    methods(Static)
        function [quadWeights, quadPoints, numIPs] = IntegrationRuleOfLine()
            %Integration of polynomials up to order 5 on the
            %reference segment (interval [-1,1])
            quadPoints = zeros(3,1);
            %xi-values (only one coordinate in reference system)
            quadPoints(1) = -sqrt(3/5);
            quadPoints(2) = 0;
            quadPoints(3) = sqrt(3/5);
            %weights
            quadWeights = zeros(3,1);
            quadWeights(1) = 5/9;
            quadWeights(2) = 8/9;
            quadWeights(3) = 5/9;
            %
            numIPs = 3;
        end
        
        function [quadWeights, quadPoints, numIPs] = IntegrationRuleOfTriangle()
            %Integration of polynomials up to order 5 on the unit
            %triangle (taken from D. Braess, Finite Elemente, Springer Verlag)
            quadPoints = zeros(7,2);
            %xi-values
            quadPoints(1,1) = 1/3;
            quadPoints(2,1) = (6 + sqrt(15))/21;
            quadPoints(3,1) = (9 - 2*sqrt(15))/21;
            quadPoints(4,1) = (6 + sqrt(15))/21;
            quadPoints(5,1) = (6 - sqrt(15))/21;
            quadPoints(6,1) = (9 + 2*sqrt(15))/21;
            quadPoints(7,1) = (6 - sqrt(15))/21;
            %eta-values
            quadPoints(1,2) = 1/3;
            quadPoints(2,2) = (6 + sqrt(15))/21;
            quadPoints(3,2) = (6 + sqrt(15))/21;
            quadPoints(4,2) = (9 - 2*sqrt(15))/21;
            quadPoints(5,2) = (6 - sqrt(15))/21;
            quadPoints(6,2) = (6 - sqrt(15))/21;
            quadPoints(7,2) = (9 + 2*sqrt(15))/21;
            %weights
            quadWeights = zeros(7,1);
            quadWeights(1) = 9/80;
            quadWeights(2:4) = (155 + sqrt(15))/2400;
            quadWeights(5:7) = (155 - sqrt(15))/2400;
            %
            numIPs = 7;
        end
        
        function u = solve(A, f)
            leaveOutDofs = find(sum(abs(A)) == 0);
            if (numel(leaveOutDofs)==0)
                u = A\f;
            else
                %first elimination of possible zero rows and columns
                A = A(any(A), :);
                A = A(:, any(A));
                %adjusting r.h.s
                f(leaveOutDofs) = [];
                %then solve linear system of equations
                
                u = A\f;
                %expand solution vector to original size and fill with zeros
                %reduce to true vector
                u = u(:);
                %do insertion
                tempPos = false(1,numel(u)+numel(leaveOutDofs));
                tempPos(leaveOutDofs) = true;
                tempVec = double(tempPos);
                tempVec(tempPos) = 0;
                tempVec(~tempPos) = u;
                u = tempVec;
            end
        end
        
        function mesh = load_gmsh(filename)
            % Reads a mesh in msh format, version 1 or 2
            
            % Copyright (C) 10/2007 R Lorph?vre (r.lorphevre@ulg.ac.be)
            
            % Based on load_gmsh supplied with gmsh-2.0 and load_gmsh2 from JP
            % Moitinho de Almeida (moitinho@civil.ist.utl.pt)
            
            % number of nodes in function of the element type
            msh.NODES_PER_TYPE_OF_ELEMENT = [
                2  3  4  4  8  6  5  3  6  9 10 27 18 14  1  8 20 15 13 ];
            
            % The format 2 don't sort the elements by reg phys but by
            % reg-elm. If this classification is important for your program,
            % use this (after calling this function):
            %
            % [OldRowNumber, NewRowNumber] = sort(OldMatrix(:,SortColumn));
            % NewMatrix = OldMatrix(NewRowNumber,:);
            %
            % Change the name of OldMatrix and NewMatrix with the name of yours
            % SortColumn by the number of the last column
            
            mesh = [];
            mesh.MIN = zeros(3, 1);
            mesh.MAX = zeros(3, 1);
            fid = fopen(filename, 'r');
            disp (' ')
            while 1
                endoffile = 0;
                while 1
                    tline = fgetl(fid);
                    if feof(fid), endoffile=1; break, end
                    if (tline(1) == '$' )
                        if tline(2) == 'N' && tline(3) == 'O'
                            fileformat = 1 ;
                            break
                        end
                        if tline(2) == 'M' && tline(3) == 'e'
                            fileformat = 2;
                            tline = fgetl(fid);
                            %disp('Mesh Type')
                            %disp (tline)
                            tline = fgetl(fid);
                            if (tline(1) == '$' && tline(2) == 'E'&& tline(3) == 'n')
                                tline = fgetl(fid);
                                break
                            else
                                disp (' This program can only read ASCII mesh file')
                                disp (' of format 1 or 2 from GMSH, try again?')
                                endoffile=1;
                                break
                            end
                        end
                        if tline(2) == 'E' && (tline(3) == 'L' || tline(3) == 'l' )
                            break
                        end
                    end
                end
                if endoffile == 1, break, end
                if tline(2) == 'N' && ( tline(3) == 'O' || tline(3) == 'o' )
                    %disp('reading nodes')
                    mesh.nbNod = fscanf(fid, '%d', 1);
                    mesh.POS = zeros(mesh.nbNod, 3);
                    for I=1:mesh.nbNod
                        iNod = fscanf(fid, '%d', 1);
                        X = fscanf(fid, '%g', 3);
                        IDS(iNod) = I;
                        if (I == 1)
                            mesh.MIN = X;
                            mesh.MAX = X;
                        else
                            if(mesh.MAX(1) < X(1)) mesh.MAX(1) = X(1); end
                            if(mesh.MAX(2) < X(2)) mesh.MAX(2) = X(2); end
                            if(mesh.MAX(3) < X(3)) mesh.MAX(3) = X(3); end
                            if(mesh.MIN(1) > X(1)) mesh.MIN(1) = X(1); end
                            if(mesh.MIN(2) > X(2)) mesh.MIN(2) = X(2); end
                            if(mesh.MIN(3) > X(3)) mesh.MIN(3) = X(3); end
                        end
                        mesh.POS(I, 1) = X(1);
                        mesh.POS(I, 2) = X(2);
                        mesh.POS(I, 3) = X(3);
                    end
                    tline = fgetl(fid);
                    %disp('nodes have been read')
                elseif tline(2) == 'E' && ( tline(3) == 'L' || tline(3) == 'l')
                    %disp('reading elements')
                    mesh.nbElm = fscanf(fid, '%d', 1);
                    if (fileformat == 1)
                        nbinfo = 5;
                        tags = 3;
                    end
                    if (fileformat == 2)
                        nbinfo = 4;
                        tags = 4;
                    end
                    mesh.ELE_INFOS = zeros(mesh.nbElm, nbinfo);
                    mesh.nbPoints = 0;
                    mesh.nbLines = 0;
                    mesh.nbTriangles = 0;
                    mesh.nbQuads = 0;
                    mesh.nbTets = 0;
                    mesh.nbHexas = 0;
                    mesh.nbPrisms = 0;
                    mesh.nbPyramids = 0;
                    %own addition%
                    %%%%%%%%%%%%%%
                    mesh.nbLines3 = 0;
                    mesh.nbTriangles6 = 0;
                    %%%%%%%%%%%%%
                    % comment next 8 lines to get "tight" arrays (will slow down reading)
                    mesh.POINTS = zeros(mesh.nbElm, 2);
                    mesh.LINES = zeros(mesh.nbElm, 3);
                    mesh.TRIANGLES = zeros(mesh.nbElm, 4);
                    mesh.QUADS = zeros(mesh.nbElm, 5);
                    mesh.TETS = zeros(mesh.nbElm, 5);
                    mesh.HEXAS = zeros(mesh.nbElm, 9);
                    mesh.PRISMS = zeros(mesh.nbElm, 7);
                    mesh.PYRAMIDS = zeros(mesh.nbElm, 6);
                    %own addition%
                    %%%%%%%%%%%%%%
                    mesh.LINES3 = zeros(mesh.nbElm, 4);
                    mesh.TRIANGLES6 = zeros(mesh.nbElm, 7);
                    %%%%%%%%%%%%%
                    for(I = 1:mesh.nbElm)
                        mesh.ELE_INFOS(I, :) = fscanf(fid, '%d', nbinfo);
                        if (fileformat == 1)
                            % take the mesh.ELE_INFOS(I, 5) nodes of the element
                            NODES_ELEM = fscanf(fid, '%d', mesh.ELE_INFOS(I, 5));
                        end
                        if (fileformat == 2)
                            mesh.ELE_TAGS(I,:) = fscanf(fid, '%d', (mesh.ELE_INFOS(I,3)-1));
                            % take the msh.NODES_PER_TYPE_OF_ELEMENT (mesh.ELE_INFOS(I, 2)) nodes of the element
                            NODES_ELEM = fscanf(fid, '%d', msh.NODES_PER_TYPE_OF_ELEMENT (mesh.ELE_INFOS(I, 2)) );
                        end
                        if(mesh.ELE_INFOS(I, 2) == 15) %% point
                            mesh.nbPoints = mesh.nbPoints + 1;
                            mesh.POINTS(mesh.nbPoints, 1) = IDS(NODES_ELEM(1));
                            mesh.POINTS(mesh.nbPoints, 2) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 1) %% line
                            mesh.nbLines = mesh.nbLines + 1;
                            mesh.LINES(mesh.nbLines, 1) = IDS(NODES_ELEM(1));
                            mesh.LINES(mesh.nbLines, 2) = IDS(NODES_ELEM(2));
                            mesh.LINES(mesh.nbLines, 3) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 2) %% triangle
                            mesh.nbTriangles = mesh.nbTriangles + 1;
                            mesh.TRIANGLES(mesh.nbTriangles, 1) = IDS(NODES_ELEM(1));
                            mesh.TRIANGLES(mesh.nbTriangles, 2) = IDS(NODES_ELEM(2));
                            mesh.TRIANGLES(mesh.nbTriangles, 3) = IDS(NODES_ELEM(3));
                            mesh.TRIANGLES(mesh.nbTriangles, 4) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 3) %% quadrangle
                            mesh.nbQuads = mesh.nbQuads + 1;
                            mesh.QUADS(mesh.nbQuads, 1) = IDS(NODES_ELEM(1));
                            mesh.QUADS(mesh.nbQuads, 2) = IDS(NODES_ELEM(2));
                            mesh.QUADS(mesh.nbQuads, 3) = IDS(NODES_ELEM(3));
                            mesh.QUADS(mesh.nbQuads, 4) = IDS(NODES_ELEM(4));
                            mesh.QUADS(mesh.nbQuads, 5) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 4) %% tetrahedron
                            mesh.nbTets = mesh.nbTets + 1;
                            mesh.TETS(mesh.nbTets, 1) = IDS(NODES_ELEM(1));
                            mesh.TETS(mesh.nbTets, 2) = IDS(NODES_ELEM(2));
                            mesh.TETS(mesh.nbTets, 3) = IDS(NODES_ELEM(3));
                            mesh.TETS(mesh.nbTets, 4) = IDS(NODES_ELEM(4));
                            mesh.TETS(mesh.nbTets, 5) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 5) %% hexahedron
                            mesh.nbHexas = mesh.nbHexas + 1;
                            mesh.HEXAS(mesh.nbHexas, 1) = IDS(NODES_ELEM(1));
                            mesh.HEXAS(mesh.nbHexas, 2) = IDS(NODES_ELEM(2));
                            mesh.HEXAS(mesh.nbHexas, 3) = IDS(NODES_ELEM(3));
                            mesh.HEXAS(mesh.nbHexas, 4) = IDS(NODES_ELEM(4));
                            mesh.HEXAS(mesh.nbHexas, 5) = IDS(NODES_ELEM(5));
                            mesh.HEXAS(mesh.nbHexas, 6) = IDS(NODES_ELEM(6));
                            mesh.HEXAS(mesh.nbHexas, 7) = IDS(NODES_ELEM(7));
                            mesh.HEXAS(mesh.nbHexas, 8) = IDS(NODES_ELEM(8));
                            mesh.HEXAS(mesh.nbHexas, 9) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 6) %% prism
                            mesh.nbPrisms = mesh.nbPrisms + 1;
                            mesh.PRISMS(mesh.nbPrisms, 1) = IDS(NODES_ELEM(1));
                            mesh.PRISMS(mesh.nbPrisms, 2) = IDS(NODES_ELEM(2));
                            mesh.PRISMS(mesh.nbPrisms, 3) = IDS(NODES_ELEM(3));
                            mesh.PRISMS(mesh.nbPrisms, 4) = IDS(NODES_ELEM(4));
                            mesh.PRISMS(mesh.nbPrisms, 5) = IDS(NODES_ELEM(5));
                            mesh.PRISMS(mesh.nbPrisms, 6) = IDS(NODES_ELEM(6));
                            mesh.PRISMS(mesh.nbPrisms, 7) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 7) %% pyramid
                            mesh.nbPyramids = mesh.nbPyramids + 1;
                            mesh.PYRAMIDS(mesh.nbPyramids, 1) = IDS(NODES_ELEM(1));
                            mesh.PYRAMIDS(mesh.nbPyramids, 2) = IDS(NODES_ELEM(2));
                            mesh.PYRAMIDS(mesh.nbPyramids, 3) = IDS(NODES_ELEM(3));
                            mesh.PYRAMIDS(mesh.nbPyramids, 4) = IDS(NODES_ELEM(4));
                            mesh.PYRAMIDS(mesh.nbPyramids, 5) = IDS(NODES_ELEM(5));
                            mesh.PYRAMIDS(mesh.nbPyramids, 6) = IDS(NODES_ELEM(6));
                            mesh.PYRAMIDS(mesh.nbPyramids, 7) = mesh.ELE_INFOS(I, tags);
                        end
                        %own addition%
                        %%%%%%%%%%%%%%
                        if(mesh.ELE_INFOS(I, 2) == 8) %% second order line
                            mesh.nbLines3 = mesh.nbLines3 + 1;
                            mesh.LINES3(mesh.nbLines3, 1) = IDS(NODES_ELEM(1));
                            mesh.LINES3(mesh.nbLines3, 2) = IDS(NODES_ELEM(2));
                            mesh.LINES3(mesh.nbLines3, 3) = IDS(NODES_ELEM(3));
                            mesh.LINES3(mesh.nbLines3, 4) = mesh.ELE_INFOS(I, tags);
                        end
                        if(mesh.ELE_INFOS(I, 2) == 9) %% second order triangle
                            mesh.nbTriangles6 = mesh.nbTriangles6 + 1;
                            mesh.TRIANGLES6(mesh.nbTriangles6, 1) = IDS(NODES_ELEM(1));
                            mesh.TRIANGLES6(mesh.nbTriangles6, 2) = IDS(NODES_ELEM(2));
                            mesh.TRIANGLES6(mesh.nbTriangles6, 3) = IDS(NODES_ELEM(3));
                            mesh.TRIANGLES6(mesh.nbTriangles6, 4) = IDS(NODES_ELEM(4));
                            mesh.TRIANGLES6(mesh.nbTriangles6, 5) = IDS(NODES_ELEM(5));
                            mesh.TRIANGLES6(mesh.nbTriangles6, 6) = IDS(NODES_ELEM(6));
                            mesh.TRIANGLES6(mesh.nbTriangles6, 7) = mesh.ELE_INFOS(I, tags);
                        end
                        %%%%%%%%%%%%%
                    end
                    tline = fgetl(fid);
                    disp('Mesh file: elements have been read')
                    
                end
            end
            
            fclose(fid);
        end
    end
end
