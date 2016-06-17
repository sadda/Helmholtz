function [x, lambda, iteration, resGibbs, mu] = SolveSemismoothEquality(A, b, C, d, x, lambda, TriInfo, matrices, lowerBound)
    if nargin < 9
        lowerBound = zeros(size(x));
    end
    
    iterMax   = 100;
    c         = 1;
    iteration = 0;
    resGibbs  = Inf;
    tol       = 1e-10;
    
    AC        = [A C'; C sparse(size(C,1),size(C,1))];
    active    = (x+c*lambda<=lowerBound);
    inactive  = ~active;
    
    while iteration <= iterMax && resGibbs >= tol
        iteration = iteration + 1;
        x         = zeros(size(x));
        lambda    = zeros(size(x));
        
        inactivePlus = logical([inactive; ones(size(C,1),1)]);
        phiMu        = AC(inactivePlus, inactivePlus) \ [-b(inactive); d];
        x(inactive)  = phiMu(1:sum(inactive));
        mu           = phiMu(sum(inactive)+1:end);
        
        if isempty(C)
            lambda(active) = -A(active,inactive)*x(inactive) - b(active);
        else
            lambda(active) = -A(active,inactive)*x(inactive) - C(:,active)'*mu - b(active);
        end
        
        active    = x+c*lambda <= lowerBound;
        inactive  = ~active;
        lambdaOpt = x - max(x + c*lambda,lowerBound);
        resGibbs  = sqrt(lambdaOpt'*lambdaOpt);
    end
    if iteration >= iterMax
        error('Careful, projection got stuck.');
    end
end