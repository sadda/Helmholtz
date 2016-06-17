function [phiProj,lambdaProj,resGibbs,iteration] = projection2Gibbs(phiCenter,phi0,matrices,lambda0,TriInfo)
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
    
    [phiProj, lambdaProj, iteration, resGibbs] = SolveSemismoothEquality(A, b, C, d, phiProj, lambdaProj, TriInfo, matrices);
    
    phiProj=reshape(phiProj,[],TriInfo.sizePhi);
end