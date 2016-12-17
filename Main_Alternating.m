% clear all;
% close all;

addpath('./OldCodes');
addpath(genpath('./P1AFEM'));

refineMesh   = 4;
drawResults  = 0;
iterMax      = 500;
iterMaxIn    = 50;
alpha        = 8e-5;



% alpha = 1e-6;


method       = 0;
refineCount  = 2;
coarsenCount = 2;

tryNumber     = 19;
regThetaL1All = linspace(1e-1, 1e-0, tryNumber);
regPhiL1All   = linspace(0, 0, tryNumber);
regPhiL2All   = linspace(0, 0, tryNumber);
cutThreshold  = 0.2;

% for tryIndex = 1:tryNumber
% for tryIndex = 1:1
    % for tryIndex = [1 5 11 19]
    % for tryIndex = [5 11 19]
    for tryIndex = 15
    tInit1      = 1e3;
    tInit2      = 1e3;
    
    if method == 1
        dirNameBase = sprintf('Results_Alternating3_RegThetaL1=%1.3f', regThetaL1All(tryIndex));
    elseif method == 0
        dirNameBase = sprintf('Results_Alternating3_Cutoff=%1.1f', cutThreshold);
    end
    picName     = fullfile(dirNameBase, 'Pictures');
    
    if exist(picName, 'dir')
        rmdir(picName, 's');
    end
    mkdir(picName);
    
    for meshIndex = 1:refineMesh
        if meshIndex == 1
            data           = load('MeshesCreated/Mesh_GeFree_AirFree/Data.mat');
            % data           = load('MeshesCreated/Mesh_AllFree/Data.mat');
            TriInfo        = data.TriInfo;
            matrices       = data.matrices;
            Transformation = data.Transformation;
            npoint0        = TriInfo.npoint;
            meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))];
            epsilon        = 2*max(meshMaxElement);
            phi            = 1/TriInfo.sizePhi*ones(TriInfo.npointRed, TriInfo.sizePhi);
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
            
            epsilon        = epsilon / 2;
            meshMaxElement = meshMaxElement / 2;
        end
        
        [constants, material]  = ObtainData(epsilon, alpha);
        material.epsilonR      = [4.2 3.4 3.4 2 1.5 1].^2;
        material.epsilonR(2:3) = [];
        
        
        
        
        
        material.epsilonR = (2*pi/1.64)^2*material.epsilonR;
        
   
        
        


        
        
        constants.regThetaL1   = regThetaL1All(tryIndex);
        constants.regPhiL1     = regPhiL1All(tryIndex);
        constants.regPhiL2     = regPhiL2All(tryIndex);
        
        fixGePrev              = fixGe;
        
        for iterIn=1:iterMaxIn
            if sum(fixGe) == 0
                error('Nono, you have to presribe the optical cavity somewhere');
            end
            
            dirName1                 = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
            dirName2                 = sprintf('Results_Ref%d_%d_Helm', meshIndex-1, iterIn);
            dirName1                 = fullfile(dirNameBase, dirName1);
            dirName2                 = fullfile(dirNameBase, dirName2);
            
            [TriInfo, matrices, phi] = ModifyMatrices(fixGe, TriInfo, matrices, Transformation, phi);
            totalArea                = TriInfo.totalArea;
            [phi, tInit1, data]      = ProjectedGradients_Old(TriInfo, Transformation, matrices, material, constants, dirName1, iterMax, drawResults, phi, tInit1);
            
            SymmetryError(phi, TriInfo, 1)
            phi = SymmetryCompute(phi, TriInfo, 1);
            
            [TriInfo, matrices, phi] = ModifyMatrices(false(TriInfo.npoint, 1), TriInfo, matrices, Transformation, phi);
            
            if method == 2
                % TODO uplne prepsat podle ostatnich.
                %             phiGe         = ProlongPhi(phi(:,1), TriInfo, 1);
                %             phiGeNorm     = sqrt(phiGe'*matrices.Mloc(1:TriInfo.npoint,1:TriInfo.npoint)*phiGe);
                %             TriInfo.phiGe = phiGe / phiGeNorm;
                %             phi2          = ProjectedGradients3(TriInfo, Transformation, matrices, material, constants, dirName, iterMax, drawResults, phi, tInit);
                
                
                wenfowefoiwef
                
            elseif method == 1
                phiProlonged                     = ProlongPhi(phi, TriInfo);
                phiGeNorm                        = sqrt(phiProlonged(:,1)'*matrices.Mloc(1:TriInfo.npoint,1:TriInfo.npoint)*phiProlonged(:,1));
                TriInfo.phiGe                    = phiProlonged(:,1) / phiGeNorm;
                [phi2, tInit2, Theta, dataEigen] = ProjectedGradients2(TriInfo, Transformation, matrices, material, constants, dirName2, iterMax, drawResults, phi, tInit2, dataEigen);
                
                SymmetryError(phi2, TriInfo, 1)
                phi2 = SymmetryCompute(phi2, TriInfo, 1);
                
                phi2Prolonged                    = ProlongPhi(phi2, TriInfo);
                fixGe                            = phi2Prolonged(:,1) > 0.99;
            else
                TriInfo.phiGe = 0;
                [~,~,~,~,~,~,~,~,~,Theta] = ComputeData2(phi, TriInfo, Transformation, matrices, constants, material, 0, []);
                fixGe = Theta >= cutThreshold*max(Theta(:));
            end
            
            sum(fixGe)
            sum(abs(fixGe - fixGePrev))
            
            %% Draw results
            
            phiProlonged = ProlongPhi(phi, TriInfo);
            
            minX      = min(TriInfo.x);
            maxX      = max(TriInfo.x);
            minY      = min(TriInfo.y);
            maxY      = max(TriInfo.y);
            
            fig = figure('visible', 'off');
            hold on;
            trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiProlonged*(1:TriInfo.sizePhi)');
            plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
            title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', totalArea, data.J1, data.J1 / totalArea));
            view(2);
            shading interp;
            colormap(colormap(gray));
            set(gca,'xcolor',get(gcf,'color'));
            set(gca,'xtick',[]);
            set(gca,'ycolor',get(gcf,'color'));
            set(gca,'ytick',[]);
            saveas(fig, fullfile(picName, sprintf('Phi_Ref%d_%d.jpg', meshIndex-1, iterIn)));
            
            fig = figure('visible', 'off');
            trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, Theta);
            title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', totalArea, data.J1, data.J1 / totalArea));
            view(2);
            shading interp;
            colormap jet;
            colorbar;
            hold on;
            set(gca,'xcolor',get(gcf,'color'));
            set(gca,'xtick',[]);
            set(gca,'ycolor',get(gcf,'color'));
            set(gca,'ytick',[]);
            % Draw boundary edges for fixGePrev
            xEle     = TriInfo.x(TriInfo.e2p);
            yEle     = TriInfo.y(TriInfo.e2p);
            fixGeEle = sum(fixGePrev(TriInfo.e2p), 2) == 3;
            xDraw = xEle(fixGeEle,:);
            xDraw = xDraw(:,[1 1 2 2 3 3]);
            xDraw = reshape(xDraw, [], 2);
            yDraw = yEle(fixGeEle,:);
            yDraw = yDraw(:,[1 1 2 2 3 3]);
            yDraw = reshape(yDraw, [], 2);
            sortIndex          = (xDraw(:,1) > xDraw(:,2)) | (xDraw(:,1) == xDraw(:,2) & yDraw(:,1) > yDraw(:,2));
            xDraw(sortIndex,:) = xDraw(sortIndex,[2 1]);
            yDraw(sortIndex,:) = yDraw(sortIndex,[2 1]);
            xyDraw   = sortrows([xDraw yDraw]);
            xyDiff1  = diff(xyDraw);
            xyDiff2  = [ones(1, 4); xyDiff1];
            xyDiff3  = [xyDiff1; ones(1, 4)];
            xyUnique = sum(abs(xyDiff2),2)~=0 & sum(abs(xyDiff3),2)~=0;
            if ~isempty(xyDraw)
                xyDraw   = xyDraw(xyUnique,:);
                plot3(xyDraw(:,1:2)', xyDraw(:,3:4)', 2*max(Theta)*ones(size(xyDraw,1),2)', 'k');
            else
                scatter3(TriInfo.x(fixGe), TriInfo.y(fixGe), 2*max(Theta)*ones(sum(fixGe),1), 'k');
            end
            % Draw boundary edges for fixGe
            xEle     = TriInfo.x(TriInfo.e2p);
            yEle     = TriInfo.y(TriInfo.e2p);
            fixGeEle = sum(fixGe(TriInfo.e2p), 2) == 3;
            xDraw = xEle(fixGeEle,:);
            xDraw = xDraw(:,[1 1 2 2 3 3]);
            xDraw = reshape(xDraw, [], 2);
            yDraw = yEle(fixGeEle,:);
            yDraw = yDraw(:,[1 1 2 2 3 3]);
            yDraw = reshape(yDraw, [], 2);
            sortIndex          = (xDraw(:,1) > xDraw(:,2)) | (xDraw(:,1) == xDraw(:,2) & yDraw(:,1) > yDraw(:,2));
            xDraw(sortIndex,:) = xDraw(sortIndex,[2 1]);
            yDraw(sortIndex,:) = yDraw(sortIndex,[2 1]);
            xyDraw   = sortrows([xDraw yDraw]);
            xyDiff1  = diff(xyDraw);
            xyDiff2  = [ones(1, 4); xyDiff1];
            xyDiff3  = [xyDiff1; ones(1, 4)];
            xyUnique = sum(abs(xyDiff2),2)~=0 & sum(abs(xyDiff3),2)~=0;
            if ~isempty(xyDraw)
                xyDraw   = xyDraw(xyUnique,:);
                plot3(xyDraw(:,1:2)', xyDraw(:,3:4)', 2*max(Theta)*ones(size(xyDraw,1),2)', 'k--');
            else
                scatter3(TriInfo.x(fixGe), TriInfo.y(fixGe), 2*max(Theta)*ones(sum(fixGe),1), 'k', 'filled');
            end
            saveas(fig, fullfile(picName, sprintf('Theta_Ref%d_%d.jpg', meshIndex-1, iterIn)));
            
            %% Stopping criterion
            
            if isequal(fixGe, fixGePrev)
                break;
            end
            fixGePrev = fixGe;
        end
    end
end


% exit;





