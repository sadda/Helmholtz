function [phiProj, lambdaProj, resGibbs, iteration] = ProjectionGibbs(phiCenter, phi0, matrices, lambda0, TriInfo)
    % Computes the H1 projection onto the Gibbs simplex via the Semismooth Newton method
    
    %% Initialize data
    
    phiProj = phi0(:);
    if isempty(lambda0)
        lambdaProj = rand(size(phiCenter(:)));
    else
        lambdaProj = lambda0(:);
    end
    
    iterMax   = 100;
    kappa     = 1;
    iteration = 0;
    resGibbs  = Inf;
    tol       = 1e-10;
    
    A         = matrices.H1scalProlong;
    b         = -matrices.H1scalProlong*phiCenter(:);
    C         = repmat(speye(sum(TriInfo.phiRowsFree)), 1, TriInfo.sizePhi);
    d         = ones(sum(TriInfo.phiRowsFree), 1);
    
    AC        = [A C'; C sparse(size(C,1),size(C,1))];
    active    = phiProj+kappa*lambdaProj <= 0;
    inactive  = ~active;
    
    %% Run Semismooth Newton method
    
    while iteration <= iterMax && resGibbs >= tol
        iteration         = iteration + 1;
        phiProj           = zeros(size(phiProj));
        lambdaProj        = zeros(size(phiProj));
        
        inactivePlus      = logical([inactive; ones(size(C,1),1)]);
        phiMu             = AC(inactivePlus, inactivePlus) \ [-b(inactive); d];
        phiProj(inactive) = phiMu(1:sum(inactive));
        mu                = phiMu(sum(inactive)+1:end);
        
        if isempty(C)
            lambdaProj(active) = -A(active,inactive)*phiProj(inactive) - b(active);
        else
            lambdaProj(active) = -A(active,inactive)*phiProj(inactive) - C(:,active)'*mu - b(active);
        end
        
        active    = phiProj+kappa*lambdaProj <= 0;
        inactive  = ~active;
        lambdaOpt = phiProj - max(phiProj + kappa*lambdaProj, 0);
        resGibbs  = sqrt(lambdaOpt'*lambdaOpt);
    end
    if iteration >= iterMax
        error('Careful, projection got stuck.');
    end
    
    phiProj = reshape(phiProj, [], TriInfo.sizePhi);
end


