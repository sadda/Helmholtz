function phiNew = ProlongPhiMesh(phiProlonged, TriInfo1, TriInfo2, normalize)
    if nargin < 4
        normalize = 1;
    end
    % Convert u from this mesh
    x1 = TriInfo1.x;
    y1 = TriInfo1.y;
    % To this mesh
    x2 = TriInfo2.x;
    y2 = TriInfo2.y;
    
    phiNew = zeros(length(x2), size(phiProlonged,2));
    for i=1:size(phiProlonged,2)
        phiNew(:,i) = griddata(x1, y1, phiProlonged(:,i), x2, y2, 'linear');
    end
    if normalize
        phiNew = max(phiNew, 0);
        phiNew = phiNew ./ repmat(sum(phiNew,2), 1, size(phiNew,2));
    end
end