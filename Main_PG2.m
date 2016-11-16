clear all;
close all;

addpath('./p1afem');
addpath(genpath('./knedlsepp-p1afem-nD-26232ef'));

tInit        = 1e3;
refineMesh   = 5;
drawResults  = 1;
drawResultsH = 1;
IterMax      = 3000;
alpha        = 1e-3;
coarsenType  = 2;
coarsenCount = 3;

for meshIndex = 1:refineMesh
    % Construct mesh and determine the starting point
    if meshIndex == 1
        load('MeshesCreated/Mesh_GeFree/Data.mat');
        meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))];
        epsilon        = 2*max(meshMaxElement);
        phi            = 1/TriInfo.sizePhi*ones(TriInfo.npointRed, TriInfo.sizePhi);
        npoint0        = TriInfo.npoint;
    else
        TriInfoOld   = TriInfo;
        phiProlonged = ProlongPhi(phi, TriInfoOld);
        
        mesh.elements    = TriInfoOld.e2p;
        mesh.coordinates = [TriInfoOld.x TriInfoOld.y];
        dirichlet        = getBoundary(mesh);
        neumann          = [];
        
        if drawResultsH
            fig = figure();
            for k=1:TriInfoOld.nelement
                patch(TriInfoOld.x(TriInfoOld.e2p(k,:))',TriInfoOld.y(TriInfoOld.e2p(k,:))',1);
                hold on;
            end
            saveas(fig, ['Mesh', int2str(meshIndex), '_1'], 'jpg');
        end
        
        % Perform refinement
        refinementTol   = 1e-6;
        markedTriangles = zeros(TriInfoOld.nelement, 1);
        for k = 1:TriInfoOld.nelement
            values = phiProlonged(TriInfoOld.e2p(k,:),:);
            if any(values(:) < 1 - refinementTol & values(:) > refinementTol) || norm(values(1,:) - values(2,:)) > refinementTol || norm(values(1,:) - values(3,:)) > refinementTol
                markedTriangles(k) = 1;
            end
        end
        [coordinates, elements, phi, dirichlet, neumann] = refineNVBModified([TriInfoOld.x TriInfoOld.y], TriInfoOld.e2p, phiProlonged, dirichlet, neumann, logical(markedTriangles));
        
        if drawResultsH
            fig = figure();
            for k=1:size(elements,1)
                patch(coordinates(elements(k,:),1)',coordinates(elements(k,:),2)',1);
                hold on;
            end
            saveas(fig, ['Mesh', int2str(meshIndex), '_2'], 'jpg');
        end
        
        coordinatesNumber = size(coordinates,1);
        switch coarsenType
            case 1
                for i = 1:coarsenCount
                    markedTriangles = zeros(size(elements,1), 1);
                    for k = 1:size(elements,1)
                        values = phi(elements(k,:),:);
                        if all(values(:) >= 1 - refinementTol | values(:) <= refinementTol) && norm(values(1,:) - values(2,:)) <= refinementTol && norm(values(1,:) - values(3,:)) <= refinementTol
                            markedTriangles(k) = 1;
                        end
                    end
                    [coordinates, elements, phi, dirichlet, neumann] = coarsenNVBModified(npoint0, coordinates, elements, phi, dirichlet, neumann, logical(markedTriangles));
                end
            case 2
                for i = 1:coarsenCount
                    markedCriterion = zeros(size(elements,1), 1);
                    for k=1:size(elements,1)
                        
                        dx1 = coordinates(elements(k,2),1) - coordinates(elements(k,1),1);
                        dy1 = coordinates(elements(k,2),2) - coordinates(elements(k,1),2);
                        
                        dx2 = coordinates(elements(k,3),1) - coordinates(elements(k,1),1);
                        dy2 = coordinates(elements(k,3),2) - coordinates(elements(k,1),2);
                        
                        
                        edet  = dx1.*dy2 - dx2.*dy1;
                        dFinv = 1/edet*[dy2, -dx2; -dy1, dx1];
                        
                        dphi = dFinv'*[-1 1 0; -1 0 1];
                        
                        slocx = ([1;0]'*dphi)*edet/2;
                        slocy = ([0;1]'*dphi)*edet/2;
                        
                        slocxPhi           = 1/edet*slocx*phi(elements(k,:),:);
                        slocyPhi           = 1/edet*slocy*phi(elements(k,:),:);
                        markedCriterion(k) = max(sqrt(slocxPhi.^2 + slocyPhi.^2));
                    end
                    markedTriangles = markedCriterion <= median(markedCriterion) / 5;
                    [coordinates, elements, phi, dirichlet, neumann] = coarsenNVBModified(npoint0, coordinates, elements, phi, dirichlet, neumann, logical(markedTriangles));
                end
        end
        fprintf('The adaptive coarsening removed %d out of %d nodes.\n\n', coordinatesNumber-size(coordinates,1), coordinatesNumber);
        
        if drawResultsH
            fig = figure();
            for k=1:size(elements,1)
                patch(coordinates(elements(k,:),1)',coordinates(elements(k,:),2)',1);
                hold on;
            end
            saveas(fig, ['Mesh', int2str(meshIndex), '_3'], 'jpg');
        end
        clear fig;
        
        
        
        %         [Transformation, TriInfo, matrices, phi] = MeshCreateProlong(TriInfo, phiProlonged, 0, 1.5*meshMaxElement, phiProlonged);
        %         if extraRefinement
        %             refinedTriangles = zeros(size(e2p,1), 1);
        %             TriInfo = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),e2p,TriInfoOld.requiredX,TriInfoOld.requiredY,TriInfoOld.sizePhi);
        %             phiProlonged = ProlongPhiMesh(phi, TriInfoOld, TriInfo);
        %             phiProlonged(~TriInfo.phiRowsFree,:) = [];
        %             phiProlonged = ProlongPhi(phiProlonged, TriInfo);
        %             for k = 1:size(e2p,1)
        %                 values = phiProlonged(e2p(k,:),:);
        %                 if any(values(:) < 1 - tol & values(:) > tol) || norm(values(1,:) - values(2,:)) > tol || norm(values(1,:) - values(3,:)) > tol % Difference in phi
        %                     if max(diff(sort(coordinates(e2p(k,:),1)))) >= 0.5*meshMaxElement(1) && max(diff(sort(coordinates(e2p(k,:),2)))) >= 0.5*meshMaxElement(2) % Element is not too small
        %                         refinedTriangles(k) = 1;
        %                     end
        %                 end
        %             end
        %             [coordinates, e2p, prolongedData] = refine(coordinates, e2p, logical(refinedTriangles), prolongedData);
        %         end
        % Recompute data
        [TriInfo, Transformation] = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),elements,TriInfoOld.requiredX,TriInfoOld.requiredY,TriInfoOld.sizePhi);
        [Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
        matrices.H1scalProlong =  TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
        
        phi(~TriInfo.phiRowsFree,:) = [];
        epsilon = epsilon / 2;
        meshMaxElement = meshMaxElement / 2;
    end
    
    dirName = ['Results_Ref', int2str(meshIndex-1)];
    try
        save(['Ref', int2str(meshIndex)]);
    end
    
    [constants, material] = ObtainData(epsilon, alpha);
    [phi, tInit] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, IterMax, drawResults, phi, tInit);
end

exit;



