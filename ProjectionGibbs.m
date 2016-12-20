function [phiProj,lambdaProj,resGibbs,iteration] = ProjectionGibbs(phiCenter,phi0,matrices,lambda0,TriInfo)
    if isempty(lambda0)
        lambdaProj = rand(size(phiCenter(:)));
    else
        lambdaProj = lambda0(:);
    end
    phiProj      = phi0(:);
    phiRowsFree  = TriInfo.phiRowsFree;
    
    A = matrices.H1scalProlong;
    b = -matrices.H1scalProlong*phiCenter(:);
    C = repmat(speye(sum(phiRowsFree)), 1, TriInfo.sizePhi);
    d = ones(sum(phiRowsFree), 1);
    
    
    
    geLevel = 0.5;
    idGe    = [matrices.Id(1:TriInfo.npoint)', sparse(1, TriInfo.npoint*(TriInfo.sizePhi-1))];
    C       = [C; idGe*TriInfo.phiProlongationMatrix6];
    d       = [d; geLevel - idGe*TriInfo.phiProlongationVector(:)];
    
    
    
    [phiProj, lambdaProj, iteration, resGibbs] = SolveSemismoothEquality(A, b, C, d, phiProj, lambdaProj, TriInfo, matrices);
    
    phiProj=reshape(phiProj,[],TriInfo.sizePhi);
end