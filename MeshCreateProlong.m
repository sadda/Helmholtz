function [Transformation, TriInfo, matrices, prolongedData] = MeshCreateProlong(TriInfoOld, phi, plotFigures, meshMaxElement, dataToProlong)
    if nargin < 5
        dataToProlong = [];
    end
    addpath('./Mesh2d v24');
    
    tol = 1e-6;
    extraRefinement = 1;
    
    % Perform refinement
    refinedTriangles = zeros(TriInfoOld.nelement, 1);
    e2p              = TriInfoOld.e2p;
    for k = 1:TriInfoOld.nelement
        values = phi(e2p(k,:),:);
        if any(values(:) < 1 - tol & values(:) > tol) || norm(values(1,:) - values(2,:)) > tol || norm(values(1,:) - values(3,:)) > tol % Difference in phi
            refinedTriangles(k) = 1;
        elseif TriInfoOld.ide(k) == 1 % Where the germanium is prescribed
            refinedTriangles(k) = 1;
        end
    end
    [coordinates, e2p, prolongedData] = refine([TriInfoOld.x TriInfoOld.y],TriInfoOld.e2p, logical(refinedTriangles), dataToProlong);
    % Do extra refinement
    if extraRefinement
        refinedTriangles = zeros(size(e2p,1), 1);
        TriInfo = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),e2p,TriInfoOld.requiredX,TriInfoOld.requiredY,TriInfoOld.sizePhi,TriInfoOld.opticalCavity);
        phiProlonged = ProlongPhiMesh(phi, TriInfoOld, TriInfo);
        phiProlonged(~TriInfo.phiRowsFree,:) = [];
        phiProlonged = ProlongPhi(phiProlonged, TriInfo);
        for k = 1:size(e2p,1)
            values = phiProlonged(e2p(k,:),:);
            if any(values(:) < 1 - tol & values(:) > tol) || norm(values(1,:) - values(2,:)) > tol || norm(values(1,:) - values(3,:)) > tol % Difference in phi
                if max(diff(sort(coordinates(e2p(k,:),1)))) >= 0.5*meshMaxElement(1) && max(diff(sort(coordinates(e2p(k,:),2)))) >= 0.5*meshMaxElement(2) % Element is not too small
                    refinedTriangles(k) = 1;
                end
            end
        end
        [coordinates, e2p, prolongedData] = refine(coordinates, e2p, logical(refinedTriangles), prolongedData);
    end
    % Recompute data
    [TriInfo, Transformation] = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),e2p,TriInfoOld.requiredX,TriInfoOld.requiredY,TriInfoOld.sizePhi,TriInfoOld.opticalCavity);
    [Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
    matrices.H1scalProlong =  TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
    
    if plotFigures
        figure();
        hold off;
        for k=1:TriInfoOld.nelement
            patch(TriInfoOld.x(TriInfoOld.e2p(k,1:3))',TriInfoOld.y(TriInfoOld.e2p(k,1:3))',1);
            hold on;
        end
        MeshCreateFigures(TriInfo);
    end
end