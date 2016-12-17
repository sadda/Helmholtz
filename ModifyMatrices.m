function [TriInfo, matrices, phi] = ModifyMatrices(fixGe, TriInfo, matrices, Transformation, phi)
    
    if nargin == 5 && size(phi,1) == TriInfo.npointRed
        phi = ProlongPhi(phi, TriInfo);
    end
    
    phiPrescribed          = TriInfo.phiProlongationVector;
    phiPrescribed(:,1)     = 0;
    fixGe                  = fixGe & sum(phiPrescribed(:,2:end), 2) == 0;
    phiPrescribed(fixGe,1) = 1;
    phiPrescribed1         = phiPrescribed(:,1);

    TrD_aa    = zeros(TriInfo.nelement,2,TriInfo.nphi);
    totalArea = 0;
    
    for k = 1:TriInfo.nelement
        if all(phiPrescribed1(TriInfo.e2p(k,:)))
            edet  = Transformation{k,1};
            slocx = Transformation{k,9};
            slocy = Transformation{k,10};
            
            TrD_aa(k,1,:) = slocx;
            TrD_aa(k,2,:) = slocy;
            totalArea     = totalArea + edet/2;
        end
    end
    
    TriInfo.phiRowsFixed           = sum(phiPrescribed, 2) > 0;
    TriInfo.phiRowsFree            = ~TriInfo.phiRowsFixed;
    TriInfo.phiRowsFree6           = repmat(TriInfo.phiRowsFree, TriInfo.sizePhi, 1);
    TriInfo.phiProlongationMatrix  = speye(TriInfo.npoint, TriInfo.npoint);
    TriInfo.phiProlongationMatrix  = TriInfo.phiProlongationMatrix(:,TriInfo.phiRowsFree);
    TriInfo.phiProlongationMatrix6 = kron(eye(TriInfo.sizePhi), TriInfo.phiProlongationMatrix);
    TriInfo.phiProlongationVector  = phiPrescribed;
    TriInfo.npointRed              = sum(TriInfo.phiRowsFree);
    TriInfo.totalArea              = totalArea;
    
    matrices.TrD           = sparse(TriInfo.indicesIElav(:),TriInfo.indicesJElav(:),TrD_aa(:));
    matrices.H1scalProlong = TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
    
    if nargin == 5 && ~isempty(phi)
        phi = phi(TriInfo.phiRowsFree,:);
    end
end




