function phi = ProlongPhi(phi, TriInfo, prolongGermanium)
    % Prolongs phi from the free nodes to the whole mesh
    
    if nargin < 3 || ~prolongGermanium
        sizePhi      = TriInfo.sizePhi;
        columnNumber = size(phi, 2);
        phi          = reshape(phi, [], sizePhi);
        phiFull      = zeros(TriInfo.npoint, sizePhi);
        for i=1:sizePhi
            phiFull(:,i) = TriInfo.phiProlongationMatrix*phi(:,i);
        end
        phi = reshape(phiFull + TriInfo.phiProlongationVector, [], columnNumber);
    else
        phi = TriInfo.phiProlongationMatrix*phi + TriInfo.phiProlongationVector(:,1);
    end
end