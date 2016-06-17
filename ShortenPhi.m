function phi = ShortenPhi(phi, TriInfo)
    columnNumber = size(phi, 2);
    phi          = reshape(phi, [], TriInfo.sizePhi);
    phiFull      = zeros(sum(TriInfo.phiRowsFree), TriInfo.sizePhi);
    for i=1:TriInfo.sizePhi
        phiFull(:,i) = TriInfo.phiProlongationMatrix'*phi(:,i);
    end
    phi = reshape(phiFull, [], columnNumber);
end
