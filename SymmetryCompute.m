function fSym = SymmetryCompute(f, TriInfo, isPhi)
    
    fOrig = f;
    if isPhi && size(f,1) > TriInfo.npoint
        f = reshape(f, [], TriInfo.sizePhi);
    end
    fReshaped = f;
    % Reduced phi
    if isPhi && size(f,1) == TriInfo.npointRed
        fAux = zeros(TriInfo.npoint, size(f,2));
        for i=1:size(f,2)
            fAux(:,i) = TriInfo.phiProlongationMatrix*f(:,i);
        end
        f = fAux;
    end
    % Theta
    if ~isPhi && size(f,1) == sum(~TriInfo.idp)
        fAux = zeros(TriInfo.npoint, size(f,2));
        fAux(~TriInfo.idp,:) = f;
        f = fAux;
    end
    
    fSym = f;
    for i=1:size(f,2)
        ind1    = find(TriInfo.x < 0);
        ind2    = find(TriInfo.x > 0);
        f1      = f(TriInfo.x < 0, i);
        f2      = f(TriInfo.x > 0, i);
        
        [g1,h1] = sort(f1);
        [q2,h2] = sort(f2);        
        fSym(ind1(h1),i) = 0.5*(g1+q2);
        fSym(ind2(h2),i) = 0.5*(g1+q2);
    end
    
    if isPhi && size(fReshaped,1) == TriInfo.npointRed
        fSym = fSym(TriInfo.phiRowsFree,:);
    end
    if ~isPhi && size(fReshaped,1) == sum(~TriInfo.idp)
        fSym = fSym(~TriInfo.idp,:);
    end
    if size(fOrig,1) > TriInfo.npoint
        fSym = fSym(:);
    end
end



