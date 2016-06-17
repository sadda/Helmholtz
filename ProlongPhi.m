function phi = ProlongPhi(phi, TriInfo)
    
    sizePhi      = TriInfo.sizePhi;
    columnNumber = size(phi, 2);
    phi          = reshape(phi, [], sizePhi);
    phiFull      = zeros(TriInfo.npoint, sizePhi);
    for i=1:sizePhi
        phiFull(:,i) = TriInfo.phiProlongationMatrix*phi(:,i);
    end
    phi = reshape(phiFull + TriInfo.phiProlongationVector, [], columnNumber);
end