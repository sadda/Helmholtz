clear all;
close all;

addpath('./OldCodes');
addpath(genpath('./P1AFEM'));

wavelength    = 1.64;
% wavelength    = 2*pi;
refineMesh    = 1;
drawResults   = 0;
iterMax       = 100;
iterMaxIn     = 20;
refineCount   = 2;
coarsenCount  = 3;
method        = 1;
standardStart = 1;

if method == 1
    pars            = [];
    pars.alpha      = repmat(linspace(1e-5, 1e-3, 6), 1, 6);
    % pars.regThetaL1 = kron(linspace(1e-4, 1e-2, 6), ones(1, 6));
    pars.regThetaL1 = kron(10.^(-4:1), ones(1, 6));
    pars.regPhiL1   = zeros(size(pars.alpha));
    pars.regPhiL2   = 1e-5*ones(size(pars.alpha));
elseif method == 0
    pars              = [];
    pars.alpha        = 1e-4;
    pars.cutThreshold = 0.2;
elseif method == -1
    pars                  = [];
    pars.jointObjectiveLp = repmat([1 2 8 2 8], 1, 3);
    pars.normalizationLp  = repmat([1 2 8 8 2], 1, 3);
    pars.alpha            = kron(linspace(1e-5, 1e-4, 3), ones(1, 5));
end

parfor parsIndex = 1:length(pars.alpha)
% for parsIndex = 1:length(pars.alpha)
    
    alpha  = pars.alpha(parsIndex);
    tInit1 = 1e3;
    tInit2 = 1e3;
    
    if method == 1
        dirNameBase = sprintf('Results_Alt_ThetaL1=%1.6f_PhiL1=%1.6f_PhiL2=%1.6f_Alpha=%1.6f', pars.regThetaL1(parsIndex), pars.regPhiL1(parsIndex), pars.regPhiL2(parsIndex), pars.alpha(parsIndex));
    elseif method == 0
        dirNameBase = sprintf('Results_Alt_Cut=%1.3f_Alpha=%1.6f', pars.cutThreshold(parsIndex), pars.alpha(parsIndex));
    elseif method == -1
        dirNameBase = sprintf('Results_Tog_Obj=%d_Norm=%d_Alpha=%1.6f', pars.jointObjectiveLp(parsIndex), pars.normalizationLp(parsIndex), pars.alpha(parsIndex));
    end
    picName     = fullfile(dirNameBase, 'Pictures');
    if exist(picName, 'dir')
        rmdir(picName, 's');
    end
    mkdir(picName);
    
    for meshIndex = 1:refineMesh
        % Construct mesh and determine the starting point
        if standardStart
            if meshIndex == 1
                % data = load('MeshesCreated/Mesh_GeFree/Data.mat');
                % load('MeshesCreated/Mesh_GeFree_AirFree/Data.mat');
                % load('MeshesCreated/Mesh_AllFree/Data.mat');
                data           = load('MeshesCreated/MeshRef_2.mat');
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
        else
            load('Results_Alternating3_Cutoff=0.2_First/Results_Ref4_23_Elas/DataAll.mat', 'TriInfo', 'phi');
            
            TriInfoOrig     = TriInfo;
            phiOrig         = ProlongPhi(phi, TriInfo);
            coordinatesOrig = [TriInfo.x TriInfo.y];
            
            minX1 = min(coordinatesOrig(phiOrig(:,1) > 1e-5, 1));
            minX2 = min(coordinatesOrig(phiOrig(:,2) > 1e-5, 1));
            maxY1 = max(coordinatesOrig(phiOrig(:,1) > 1e-5, 2));
            maxY2 = max(coordinatesOrig(phiOrig(:,2) > 1e-5, 2));
            
            load('MeshesCreated/Mesh_GeFree/Data.mat', 'TriInfo');
            
            mesh          = struct('coordinates', [TriInfo.x TriInfo.y], 'elements', TriInfo.e2p);
            coordinates   = mesh.coordinates;
            elements      = mesh.elements;
            dirichlet     = getBoundary(mesh);
            neumann       = [];
            overfill1     = 0.02;
            overfill2     = 0.25;
            
            % Perform refinement
            refinementTol = 1e-6;
            for i=1:(refineMesh-1)
                markedTriangles = zeros(size(elements,1), 1);
                for k = 1:size(elements,1)
                    coorX = coordinates(elements(k,:),1);
                    coorY = coordinates(elements(k,:),2);
                    if (all(coorX>=min(minX1,minX2)-overfill1) && all(coorX<=-min(minX1,minX2)+overfill1) && all(coorY>=1) && all(coorY<=max(maxY1,maxY2)+overfill1)) || (all(coorY>=1-overfill2) && all(coorY<=1+overfill2))
                        markedTriangles(k) = 1;
                    end
                end
                [coordinates, elements, ~, dirichlet, neumann] = refineNVBModified(coordinates, elements, [], dirichlet, neumann, logical(markedTriangles));
                overfill2 = overfill2 / 2;
            end
            
            % Recompute data
            [TriInfo, Transformation]           = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),elements,TriInfo.requiredX,TriInfo.requiredY,TriInfo.sizePhi);
            [Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
            matrices.H1scalProlong              = TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
            
            phi = zeros(TriInfo.npoint, TriInfo.sizePhi);
            for i=1:TriInfo.sizePhi
                phiInterp = scatteredInterpolant(TriInfoOrig.x, TriInfoOrig.y, phiOrig(:,i));
                phi(:,i)  = phiInterp(TriInfo.x, TriInfo.y);
            end
            phi = SymmetryCompute(phi, TriInfo, 1);
            phi = phi ./ repmat(sum(phi,2), 1, TriInfo.sizePhi);
            
            meshMaxElement = [min(diff(unique(TriInfo.x))); min(diff(unique(TriInfo.y)))];
            epsilon        = 2*max(meshMaxElement);
            npoint0        = TriInfo.npoint;
            fixGe          = phi(:,1) >= 0.5;
            dataEigen      = [];
            phi            = phi(sum(TriInfo.phiProlongationVector,2) == 0,:);
        end
        
        [constants, material] = ObtainData(epsilon, alpha, wavelength);
        
        if method ~= -1
            constants.regThetaL1 = pars.regThetaL1(parsIndex);
            constants.regPhiL1   = pars.regPhiL1(parsIndex);
            constants.regPhiL2   = pars.regPhiL2(parsIndex);
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
                
                phi = SymmetryCompute(phi, TriInfo, 1, 1, 1e-8);
                
                [TriInfo, matrices, phi] = ModifyMatrices(false(TriInfo.npoint, 1), TriInfo, matrices, Transformation, phi);
                
                if method == 1
                    phiProlonged                     = ProlongPhi(phi, TriInfo);
                    phiGeNorm                        = sqrt(phiProlonged(:,1)'*matrices.Mloc(1:TriInfo.npoint,1:TriInfo.npoint)*phiProlonged(:,1));
                    TriInfo.phiGe                    = phiProlonged(:,1) / phiGeNorm;
                    [phi2, tInit2, Theta, dataEigen] = ProjectedGradients2(TriInfo, Transformation, matrices, material, constants, dirName2, iterMax, drawResults, phi, tInit2, dataEigen);
                    
                    SymmetryError(phi2, TriInfo, 1)
                    phi2 = SymmetryCompute(phi2, TriInfo, 1);
                    
                    phi2Prolonged                    = ProlongPhi(phi2, TriInfo);
                    fixGe                            = phi2Prolonged(:,1) > 0.9;
                else
                    TriInfo.phiGe = 0;
                    
                    options = struct('computeG', 0, 'computeU', 1, 'symmetrize', 0, 'separateObjective', 0);
                    [~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
                    fixGe = Theta >= pars.cutThreshold(parsIndex)*max(Theta(:));
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
            
            options = [];
            options = SetField(options, 'computeU', 1);
            options = SetField(options, 'computeG', 0);
            options = SetField(options, 'symmetrize', 0);
            options = SetField(options, 'separateObjective', 0);
            options = SetField(options, 'jointObjectiveLp', pars.jointObjectiveLp(parsIndex));
            options = SetField(options, 'normalizationLp', pars.normalizationLp(parsIndex));
            
            [phi, tInit1]             = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, iterMax, drawResults, phi, tInit1, options);
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




exit;









