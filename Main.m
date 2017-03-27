clear all;
close all;

addpath(genpath('./P1AFEM'));

wavelength    = 1.64;
refineMesh    = 6;
drawResults   = 0;
iterMax       = 500;
iterMaxIn     = 50;
refineCount   = 3;
coarsenCount  = 2;
method        = 0;

%% Set parameters

if method == 2
    pars            = [];
    pars.alpha      = linspace(1e-5, 1e-3, 6);
    pars.regThetaL1 = 10.^(-4:1);
    pars.regPhiL1   = zeros(size(pars.alpha));
    pars.regPhiL2   = 1e-5*ones(size(pars.alpha));
    parsAll         = combvec(unique(pars.alpha), pars.regThetaL1, pars.regPhiL1, pars.regPhiL2);
elseif method == 1
    pars              = [];
    pars.alpha        = linspace(1e-5, 1e-4, 8);
    pars.cutThreshold = 0.2;
    parsAll           = combvec(unique(pars.alpha), pars.cutThreshold);
elseif method == 0
    pars                       = [];
    pars.alpha                 = linspace(5e-5, 5e-4, 8);
    pars.jointObjectiveThetaLp = 2;
    pars.jointObjectivePhiLp   = 1;
    pars.geLevel2              = 0;
    pars.gammaPen              = 0;
    pars.alphaGe               = 0;
    parsAll                    = combvec(unique(pars.alpha), pars.jointObjectiveThetaLp, pars.jointObjectivePhiLp, pars.geLevel2, unique(pars.gammaPen), unique(pars.alphaGe));
end

%% Run a loop across all parameter sets

% parfor parsIndex = 1:size(parsAll, 2)
for parsIndex = 1:size(parsAll, 2)
    
    options              = [];
    options.originalEps0 = 0;
    options.computeG     = 0;
    options.symmetrize   = 0;
    options.method       = method;
    
    alpha  = parsAll(1,parsIndex);
    tInit1 = 1e3;
    tInit2 = 1e3;
    
    % Copy parameters based on used method
    if method == 2
        dirNameBase = sprintf('Res_M2_%1.6f_%1.6f_%1.6f_%1.6f', alpha, parsAll(2,parsIndex), parsAll(3,parsIndex), parsAll(4,parsIndex));
    elseif method == 1
        cutThreshold = parsAll(2,parsIndex);
        dirNameBase  = sprintf('Res_M1_%1.6f_%1.3f', alpha, cutThreshold);
    elseif method == 0
        options.jointObjectiveThetaLp = parsAll(2,parsIndex);
        options.jointObjectivePhiLp   = parsAll(3,parsIndex);
        options.normalizationLp       = 2;
        options.geLevel2              = parsAll(4,parsIndex);
        options.gammaPen              = parsAll(5,parsIndex);
        options.alphaGe               = parsAll(6,parsIndex);
        dirNameBase = sprintf('Res_M0_%1.6f_%d_%d_%1.4f_%1.4f_%1.4f', alpha, options.jointObjectiveThetaLp, options.jointObjectivePhiLp, options.geLevel2, options.gammaPen, options.alphaGe);
    end
    % Create folders where the results will be saved
    picName     = fullfile(dirNameBase, 'Pictures');
    if exist(picName, 'dir')
        rmdir(picName, 's');
    end
    mkdir(picName);
    
    %% Run a loop across all meshes
    
    for meshIndex = 1:refineMesh
        
        %%  Construct initial mesh and set the starting point or refine it
        
        if meshIndex == 1
            data           = load('MeshesCreated/Mesh_New2/Data.mat');
            TriInfo        = data.TriInfo;
            Transformation = data.Transformation;
            matrices       = data.matrices;
            meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))];
            
            epsilon        = 2*max(meshMaxElement);
            phi            = 1/TriInfo.sizePhi*ones(TriInfo.npointRed, TriInfo.sizePhi);
            npoint0        = TriInfo.npoint;
            % TODO toto jeste prepis a uvidis, co to udela
            fixGe          = TriInfo.x >= -1 & TriInfo.x <= 1 & TriInfo.y >= 1.25 & TriInfo.y <= 1.5; % Used only for methods 1 and 2
            dataEigen      = [];
        else % Refine mesh
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
                    dataEigenPart = dataEigen(dataEigenKey);
                    dataToProlong = [dataToProlong dataEigenPart{2}];
                end
            end
            
            % Perform refinement refineCount times
            phi           = ProlongPhi(phi, TriInfo);
            refinementTol = 1e-6;
            for i=1:refineCount
                markedTriangles = zeros(size(elements,1), 1);
                for k = 1:size(elements,1)
                    values = phi(elements(k,:),:);
                    % Mark all triangles in the interfacial interface
                    if any(values(:) < 1 - refinementTol & values(:) > refinementTol) || norm(values(1,:) - values(2,:)) > refinementTol || norm(values(1,:) - values(3,:)) > refinementTol
                        side1 = max(diff(sort(coordinates(elements(k,:),1))));
                        side2 = max(diff(sort(coordinates(elements(k,:),2))));
                        % But not thouse who are too small
                        if side1 >= 0.9*meshMaxElement(1) && side2 >= 0.9*meshMaxElement(2)
                            markedTriangles(k) = 1;
                        end
                    end
                end
                [coordinates, elements, dataToProlong, dirichlet, neumann] = refineNVBModified(coordinates, elements, dataToProlong, dirichlet, neumann, logical(markedTriangles));
                phi = dataToProlong(:,1:TriInfo.sizePhi);
            end
            
            % Perform coarsening coarsenCount times
            coordinatesNumber = size(coordinates,1);
            for i=1:coarsenCount
                markedTriangles = zeros(size(elements,1), 1);
                for k = 1:size(elements,1)
                    values = phi(elements(k,:),:);
                    % Mark all elements outside of the interface
                    if all(values(:) >= 1 - refinementTol | values(:) <= refinementTol) && norm(values(1,:) - values(2,:)) <= refinementTol && norm(values(1,:) - values(3,:)) <= refinementTol
                        markedTriangles(k) = 1;
                    end
                end
                if any(markedTriangles)
                    [coordinates, elements, dataToProlong, dirichlet, neumann] = coarsenNVBModified(npoint0, coordinates, elements, dataToProlong, dirichlet, neumann, logical(markedTriangles));
                    phi = dataToProlong(:,1:TriInfo.sizePhi);
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
                counter = TriInfo.sizePhi+1; % Skip phi and fixGe
                for dataEigenKey = dataEigenKeys
                    dataEigenPart           = dataEigen(dataEigenKey);
                    dataEigen(dataEigenKey) = {dataEigenPart{1}, dataToProlong(:, counter+1:counter+TriInfo.sizePhi)};
                    counter                 = counter+TriInfo.sizePhi;
                end
            end
            
            epsilon = epsilon / 2;
            meshMaxElement = meshMaxElement / 2;
        end
        
        %% Obtain material properties
        
        [constants, material] = ObtainData(epsilon, alpha, wavelength, options);
        
        %% Run an optimization method
        
        if method ~= 0
            %% The alternating methods, where we look for a fixed point of the optical cavity
            
            fixGePrev          = fixGe;
            
            % Run the alternative procedure
            for iterIn=1:iterMaxIn
                if sum(fixGe) == 0
                    error('Nono, you have to presribe the optical cavity somewhere');
                end
                
                dirName1                 = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
                dirName2                 = sprintf('Results_Ref%d_%d_Helm', meshIndex-1, iterIn);
                dirName1                 = fullfile(dirNameBase, dirName1);
                dirName2                 = fullfile(dirNameBase, dirName2);
                
                % Some matrices need to be modified if the optical cavity (and thus prescribed phi) is changed
                [TriInfo, matrices, phi]          = ModifyMatrices(fixGe, TriInfo, matrices, Transformation, phi);
                options.method                    = 1;                
                [phi, tInit1, ~, dataEigen, data] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName1, iterMax, drawResults, phi, tInit1, options, dataEigen);                
                options.method                    = method;                
                data.totalArea                    = TriInfo.totalArea;
                
                phi                      = SymmetryCompute(phi, TriInfo, 1, 1, 1e-8);
                [TriInfo, matrices, phi] = ModifyMatrices(false(TriInfo.npoint, 1), TriInfo, matrices, Transformation, phi);
                
                if method == 2
                    options.regThetaL1 = parsAll(2,parsIndex);
                    options.regPhiL1   = parsAll(3,parsIndex);
                    options.regPhiL2   = parsAll(4,parsIndex);

                    phiProlonged                     = ProlongPhi(phi, TriInfo);
                    phiGeNorm                        = sqrt(phiProlonged(:,1)'*matrices.Mloc(1:TriInfo.npoint,1:TriInfo.npoint)*phiProlonged(:,1));
                    TriInfo.phiGe                    = phiProlonged(:,1) / phiGeNorm;
                    [phi2, tInit2, Theta, dataEigen] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName2, iterMax, drawResults, phi, tInit2, options, dataEigen);
                    phi2                             = SymmetryCompute(phi2, TriInfo, 1, 1, 1e-8);
                    phi2Prolonged                    = ProlongPhi(phi2, TriInfo);
                    fixGe                            = phi2Prolonged(:,1) > 0.9;
                else
                    TriInfo.phiGe             = 0;
                    [~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
                    fixGe                     = Theta >= cutThreshold*max(Theta(:));
                end
                
                fprintf('Changed nodes = %d\n', sum(abs(fixGe - fixGePrev)));
                fprintf('Total nodes   = %d\n\n', sum(abs(fixGe)));
                
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
            %% Use optimization all at once (Dirk's approach)
            
            iterIn  = 0;
            dirName = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
            dirName = fullfile(dirNameBase, dirName);
            
            [phi, tInit1, Theta, dataEigen] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, iterMax, drawResults, phi, tInit1, options, dataEigen);
            
            figure;
            phiPr = ProlongPhi(phi, TriInfo);
            trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiPr(:,1));
            view(3);
            
            figure;
            phiPr = ProlongPhi(phi, TriInfo);
            trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiPr(:,2));
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









