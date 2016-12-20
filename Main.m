% clear all;
% close all;

addpath('./OldCodes');
addpath(genpath('./P1AFEM'));


refineMesh   = 1;
drawResults  = 1;
% iterMax      = 1000;
iterMax      = 10;
iterMaxIn    = 50;
refineCount  = 2;
coarsenCount = 3;
method       = -1;

regThetaL1   = 1e-1;
regPhiL1     = 1e-5;
regPhiL2     = 1e-5;
cutThreshold = 0.2;

alphaAll = linspace(1e-4, 8e-4, 1);


for alphaIndex = 1:length(alphaAll)
    
    alpha = alphaAll(alphaIndex);
    tInit1       = 1e3;
    tInit2       = 1e3;
    
    if method == 1
        dirNameBase = sprintf('Results_Alternating3_RegThetaL1=%1.3f', regThetaL1);
    elseif method == 0
        dirNameBase = sprintf('Results_Alternating3_Cutoff=%1.1f', cutThreshold);
    elseif method == -1
        dirNameBase = sprintf('Results_Together_Alpha=%1.5f', alpha);
    end
    picName     = fullfile(dirNameBase, 'Pictures');
    if exist(picName, 'dir')
        rmdir(picName, 's');
    end
    mkdir(picName);
    
    for meshIndex = 1:refineMesh
        % Construct mesh and determine the starting point
        if meshIndex == 1
            % load('MeshesCreated/Mesh_GeFree/Data.mat');
            % load('MeshesCreated/Mesh_GeFree_AirFree/Data.mat');
            % load('MeshesCreated/Mesh_AllFree/Data.mat');
            data           = load('MeshesCreated/MeshRef_0.mat');
            TriInfo        = data.TriInfo;
            Transformation = data.Transformation;
            matrices       = data.matrices;
            meshMaxElement = [min(diff(unique(TriInfo.x))); min(diff(unique(TriInfo.y)))];
            
            epsilon        = 2*max(meshMaxElement);
            phi            = 1/TriInfo.sizePhi*ones(TriInfo.npointRed, TriInfo.sizePhi);
            npoint0        = TriInfo.npoint;
            fixGe          = TriInfo.x >= -1 & TriInfo.x <= 1 & TriInfo.y >= 1.25 & TriInfo.y <= 1.5;
            dataEigen      = [];
        else
            TriInfoOld   = TriInfo;
            phiProlonged = ProlongPhi(phi, TriInfoOld);
            
            mesh         = struct('coordinates', [TriInfoOld.x TriInfoOld.y], 'elements', TriInfoOld.e2p);
            coordinates  = mesh.coordinates;
            elements     = mesh.elements;
            dirichlet    = getBoundary(mesh);
            neumann      = [];
            
            % Collect data which need to be prolonged to the new mesh
            dataToProlong = [phiProlonged fixGe];
            if ~isempty(dataEigen)
                dataEigenKeys = keys(dataEigen);
                dataEigenKeys = [dataEigenKeys{:}];
                dataEigenKeys = dataEigenKeys(dataEigenKeys >= 0);
                for dataEigenKey = dataEigenKeys
                    data          = dataEigen(dataEigenKey);
                    id1D          = ~TriInfo.idp(1:TriInfo.npoint)==1;
                    yEigen        = zeros(TriInfo.npoint, 1);
                    yEigen(id1D)  = data{1};
                    dataToProlong = [dataToProlong yEigen data{3}];
                end
            end
            
            % Perform refinement
            phi           = ProlongPhi(phi, TriInfo);
            refinementTol = 1e-6;
            for i=1:refineCount
                markedTriangles = zeros(size(elements,1), 1);
                for k = 1:size(elements,1)
                    values = phi(elements(k,:),:);
                    % TODO rozmysli si, jestli nechces pridat automat kdyz jsme ve fixGe
                    if any(values(:) < 1 - refinementTol & values(:) > refinementTol) || norm(values(1,:) - values(2,:)) > refinementTol || norm(values(1,:) - values(3,:)) > refinementTol
                        if max(diff(sort(coordinates(elements(k,:),1)))) >= meshMaxElement(1)-refinementTol && max(diff(sort(coordinates(elements(k,:),2)))) >= meshMaxElement(2)-refinementTol % Element is not too small
                            markedTriangles(k) = 1;
                        end
                    end
                end
                [coordinates, elements, dataToProlong, dirichlet, neumann] = refineNVBModified(coordinates, elements, dataToProlong, dirichlet, neumann, logical(markedTriangles));
                phi   = dataToProlong(:,1:TriInfo.sizePhi);
            end
            
            % Perform coarsening
            coordinatesNumber = size(coordinates,1);
            for i=1:coarsenCount
                markedTriangles = zeros(size(elements,1), 1);
                for k = 1:size(elements,1)
                    values = phi(elements(k,:),:);
                    if all(values(:) >= 1 - refinementTol | values(:) <= refinementTol) && norm(values(1,:) - values(2,:)) <= refinementTol && norm(values(1,:) - values(3,:)) <= refinementTol
                        markedTriangles(k) = 1;
                    end
                end
                if any(markedTriangles)
                    [coordinates, elements, dataToProlong, dirichlet, neumann] = coarsenNVBModified(npoint0, coordinates, elements, dataToProlong, dirichlet, neumann, logical(markedTriangles));
                    phi   = dataToProlong(:,1:TriInfo.sizePhi);
                end
            end
            fprintf('The adaptive coarsening removed %d out of %d nodes.\n\n', coordinatesNumber-size(coordinates,1), coordinatesNumber);
            
            % Recompute data
            [TriInfo, Transformation]           = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),elements,TriInfoOld.requiredX,TriInfoOld.requiredY,TriInfoOld.sizePhi);
            [Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
            matrices.H1scalProlong              = TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
            
            phi   = dataToProlong(:,1:TriInfo.sizePhi);
            phi(~TriInfo.phiRowsFree,:) = [];
            fixGe = dataToProlong(:,TriInfo.sizePhi+1);
            fixGe = fixGe >= 1-1e-5; % Do not forget that you need binary values for fixGe (it is lost after prolongation)
            if ~isempty(dataEigen)
                counter = TriInfo.sizePhi+2;
                for dataEigenKey = dataEigenKeys
                    data                    = dataEigen(dataEigenKey);
                    id1D                    = ~TriInfo.idp(1:TriInfo.npoint)==1;
                    yEigen                  = dataToProlong(:,counter);
                    dataEigen(dataEigenKey) = {yEigen(id1D), data{2}, dataToProlong(:, counter+1:counter+TriInfo.sizePhi)};
                    counter                 = counter+1+TriInfo.sizePhi;
                end
            end
            
            epsilon = epsilon / 2;
            meshMaxElement = meshMaxElement / 2;
        end
        
        [constants, material] = ObtainData(epsilon, alpha, 1.64);
        
        if method ~= -1
            constants.regThetaL1 = regThetaL1;
            constants.regPhiL1   = regPhiL1;
            constants.regPhiL2   = regPhiL2;
            fixGePrev            = fixGe;
            
            for iterIn=1:iterMaxIn
                if sum(fixGe) == 0
                    error('Nono, you have to presribe the optical cavity somewhere');
                end
                
                dirName1                 = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
                dirName2                 = sprintf('Results_Ref%d_%d_Helm', meshIndex-1, iterIn);
                dirName1                 = fullfile(dirNameBase, dirName1);
                dirName2                 = fullfile(dirNameBase, dirName2);
                
                [TriInfo, matrices, phi] = ModifyMatrices(fixGe, TriInfo, matrices, Transformation, phi);
                [phi, tInit1, data]      = ProjectedGradients_Old(TriInfo, Transformation, matrices, material, constants, dirName1, iterMax, drawResults, phi, tInit1);
                data.totalArea           = TriInfo.totalArea;
                
                SymmetryError(phi, TriInfo, 1)
                phi = SymmetryCompute(phi, TriInfo, 1);
                
                [TriInfo, matrices, phi] = ModifyMatrices(false(TriInfo.npoint, 1), TriInfo, matrices, Transformation, phi);
                
                if method == 1
                    phiProlonged                     = ProlongPhi(phi, TriInfo);
                    phiGeNorm                        = sqrt(phiProlonged(:,1)'*matrices.Mloc(1:TriInfo.npoint,1:TriInfo.npoint)*phiProlonged(:,1));
                    TriInfo.phiGe                    = phiProlonged(:,1) / phiGeNorm;
                    [phi2, tInit2, Theta, dataEigen] = ProjectedGradients2(TriInfo, Transformation, matrices, material, constants, dirName2, iterMax, drawResults, phi, tInit2, dataEigen);
                    
                    SymmetryError(phi2, TriInfo, 1)
                    phi2 = SymmetryCompute(phi2, TriInfo, 1);
                    
                    phi2Prolonged                    = ProlongPhi(phi2, TriInfo);
                    fixGe                            = phi2Prolonged(:,1) > 0.3;
                else
                    TriInfo.phiGe = 0;
                    
                    options = struct('computeG', 0, 'computeU', 1, 'symmetrize', 0, 'separateObjective', 0);
                    [~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
                    fixGe = Theta >= cutThreshold*max(Theta(:));
                end
                
                sum(fixGe)
                sum(abs(fixGe - fixGePrev))
                
                
                data.iterIn    = iterIn;
                data.meshIndex = meshIndex;
                data.drawFixGe = 1;
                data.fixGe     = fixGe;
                data.fixGePrev = fixGePrev;
                
                DrawResults(phi, Theta, TriInfo, picName, data);
                
                %% Stopping criterion
                
                if isequal(fixGe, fixGePrev)
                    break;
                end
                fixGePrev = fixGe;
            end
        else
            iterIn  = 0;
            dirName = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
            dirName = fullfile(dirNameBase, dirName);
            
            [phi, tInit1] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, iterMax, drawResults, phi, tInit1);
            
            options = struct('computeG', 0, 'computeU', 1, 'symmetrize', 0, 'separateObjective', 0);
            [~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
            
            figure;
            phiPr = ProlongPhi(phi, TriInfo);
            trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiPr(:,1));
            view(3);
            
            data.J1        = NaN;
            data.totalArea = NaN;
            data.iterIn    = iterIn;
            data.meshIndex = meshIndex;
            data.drawFixGe = 0;
            
            DrawResults(phi, Theta, TriInfo, picName, data);
        end
    end
end




% exit;









