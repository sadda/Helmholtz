function phi = ShortenPhi(phi, TriInfo)
    % Multiplies phi on the full mesh by P'. Used for derivative computations
    
    columnNumber = size(phi, 2);
    phi          = reshape(phi, [], TriInfo.sizePhi);
    phiFull      = zeros(sum(TriInfo.phiRowsFree), TriInfo.sizePhi);
    for i=1:TriInfo.sizePhi
        phiFull(:,i) = TriInfo.phiProlongationMatrix'*phi(:,i);
    end
    phi = reshape(phiFull, [], columnNumber);
end
