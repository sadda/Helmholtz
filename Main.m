clear all;
close all;

addpath(genpath('./P1AFEM'));

%% Set parameters

wavelength    = 1.64;
refineMesh    = 3;                            % Number of the meshes on which the code will be run
drawResults1  = 1;                            % Draw displacement u, u_x, u_y, biaxial strain, Theta and possibly something more
drawResults2  = 1;                            % Draw only u and Theta
iterMax       = 20;                          % Maximal iteration number for the projected gradients
iterMaxIn     = 20;                           % Number of the maximal number of outer iterations for the alternating method
refineCount   = 3;                            % How many times nodes should be added under each mesh refinement
coarsenCount  = 3;                            % How many times nodes should be removed under each mesh refinement
method        = 0;                            % Used method: 0 = joint approach. 1 = alternating minimization

% Specify the parameters. If you add multiple values for one parameter, all possible combinations will be created.
if method == 1
    pars              = [];
    pars.alpha        = 7e-5;                 % Weight of the Ginzburg-Landau energy
    pars.cutThreshold = 0.2;                  % Parameter for determining the new optical cavity (higher values mean smaller optical cavity)/
    parsAll           = combvec(unique(pars.alpha), pars.cutThreshold);
elseif method == 0
    pars                       = [];
    pars.alpha                 = 4e-4;        % Weight of the Ginzburg-Landau energy
    pars.jointObjectiveThetaLp = 2;           % Power of Theta in the objective
    pars.jointObjectivePhiLp   = 1;           % Power of phi_{Ge} in the objective
    parsAll                    = combvec(unique(pars.alpha), pars.jointObjectiveThetaLp, pars.jointObjectivePhiLp);
end

%% Run a loop across all parameter sets

% parfor parsIndex = 1:size(parsAll, 2)
for parsIndex = 1:size(parsAll, 2)
    
    options              = [];
    options.originalEps0 = 0;                 % Determines whether the original elasticity should be used (where the eigenstrain is placed)
    options.computeG     = 0;                 % Determines whether the derivate should be computed. Is rewritten many times during the optimization
    options.symmetrize   = 1;                 % Performs symmetrization after computation of the gradient and Theta
    options.method       = method;
    
    alpha = parsAll(1,parsIndex);
    tInit = 1e3;                              % Initial stepsize for the projected method
    
    % Copy parameters based on used method
    if method == 1
        cutThreshold = parsAll(2,parsIndex);
        dirNameBase  = sprintf('Res_M1_%1.6f_%1.3f', alpha, cutThreshold);
    elseif method == 0
        options.jointObjectiveThetaLp = parsAll(2,parsIndex);
        options.jointObjectivePhiLp   = parsAll(3,parsIndex);
        options.normalizationLp       = options.jointObjectiveThetaLp;
        dirNameBase = sprintf('Res_M0_%1.6f_%d_%d', alpha, options.jointObjectiveThetaLp, options.jointObjectivePhiLp);
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
        
        if meshIndex == 1 % Load the initial mesh
            data           = load('MeshesCreated/Mesh_Final/Data.mat'); % Needs to be written in the way due to parfor
            TriInfo        = data.TriInfo;
            Transformation = data.Transformation;
            matrices       = data.matrices;
            meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))]; % Determines the size of the element on the interface
            
            epsilon        = 2*max(meshMaxElement); % Parameter related to the interfacial thickness. Will be decreased during the optimization
            phi            = 1/TriInfo.sizePhi*ones(TriInfo.npointRed, TriInfo.sizePhi); % Starting with uniform phi
            npoint0        = TriInfo.npoint; % Number of initial nodes. The coarsening procedure cannot remove nodes with lower number
            fixGe          = TriInfo.x >= -1 & TriInfo.x <= 1 & TriInfo.y >= 1.01 & TriInfo.y <= 1.49; % The original optical cavity for method 1
            dataEigen      = [];
        else % Refine mesh
            TriInfoOld   = TriInfo;
            phiProlonged = ProlongPhi(phi, TriInfoOld);
            
            mesh         = struct('coordinates', [TriInfoOld.x TriInfoOld.y], 'elements', TriInfoOld.e2p); % The mesh need to be put into the structure used by P1AFEM
            coordinates  = mesh.coordinates;
            elements     = mesh.elements;
            dirichlet    = getBoundary(mesh);
            neumann      = [];
            
            % Collect data which need to be prolonged to the new mesh
            dataToProlong = [phiProlonged fixGe]; % Stores the data to be prolonged in the following order: phi, phiGe, Theta_1, Theta_2, ...
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
            refinementTol = 1e-6; % For pure phases we consider inteval (refinementTol, 1-refinementTols)
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
                % Recompute the new mesh (without any matrices, they will be added later)
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
                    % Recompute the new mesh (without any matrices, they will be added later)
                    [coordinates, elements, dataToProlong, dirichlet, neumann] = coarsenNVBModified(npoint0, coordinates, elements, dataToProlong, dirichlet, neumann, logical(markedTriangles));
                    phi = dataToProlong(:,1:TriInfo.sizePhi);
                end
            end
            fprintf('The adaptive coarsening removed %d out of %d nodes.\n\n', coordinatesNumber-size(coordinates,1), coordinatesNumber);
            
            % Recompute data including matrices. Written in the simplest (read stupidest) possible way
            [TriInfo, Transformation]           = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),elements,TriInfoOld.requiredX,TriInfoOld.requiredY,TriInfoOld.sizePhi);
            [Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
            matrices.H1scalProlong              = TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
            
            % Recover data from dataToProlong
            phi   = dataToProlong(:,1:TriInfo.sizePhi);
            phi(~TriInfo.phiRowsFree,:) = []; % Phi is stored only on the reduced mesh (where no materials are prescribed)
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
            
            % Decrease the parameter epsilon by one half
            epsilon        = epsilon / 2;
            meshMaxElement = meshMaxElement / 2;
        end
        
        %% Obtain constants and material properties
        
        [constants, material] = ObtainData(epsilon, alpha, wavelength, options);
        
        %% Run an optimization method
        
        if method == 1
            %% The alternating methods, where we look for a fixed point of the optical cavity
            
            fixGePrev          = fixGe;
            
            % Run the alternative procedure
            for iterIn=1:iterMaxIn
                if sum(fixGe) == 0
                    error('Nono, you have to presribe the optical cavity somewhere');
                end
                
                dirName1 = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
                dirName2 = sprintf('Results_Ref%d_%d_Helm', meshIndex-1, iterIn);
                dirName1 = fullfile(dirNameBase, dirName1);
                dirName2 = fullfile(dirNameBase, dirName2);
                
                % Some matrices need to be modified if the optical cavity (and thus prescribed phi) is changed
                [TriInfo, matrices, phi]   = ModifyMatrices(fixGe, TriInfo, matrices, Transformation, phi);
                [phi, tInit, ~, dataEigen] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName1, iterMax, drawResults1, phi, tInit, options, dataEigen);
                options.method             = method;
                phi                        = SymmetryCompute(phi, TriInfo, 1, 1, 1e-8);
                [TriInfo, matrices, phi]   = ModifyMatrices(false(TriInfo.npoint, 1), TriInfo, matrices, Transformation, phi);
                
                % Compute Theta and the new optical cavity
                TriInfo.phiGe             = 0;
                options.method            = 0;
                [~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
                options.method            = 1;
                fixGe                     = Theta >= cutThreshold*max(Theta(:));
                
                fprintf('Changed nodes = %d\n', sum(abs(fixGe - fixGePrev)));
                fprintf('Total nodes   = %d\n\n', sum(abs(fixGe)));
                
                % Prepare data for drawing and draw it
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
            %% The joint approach
            
            iterIn  = 0;
            dirName = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
            dirName = fullfile(dirNameBase, dirName);
            
            % Run the projected gradients
            [phi, tInit, Theta, dataEigen] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, iterMax, drawResults1, phi, tInit, options, dataEigen);
            
            % And finally draw phi and Theta
            if drawResults2
                data.iterIn    = iterIn;
                data.meshIndex = meshIndex;
                data.drawFixGe = 0;
                
                DrawResults(phi, Theta, TriInfo, picName, data);
            end
        end
    end
end










