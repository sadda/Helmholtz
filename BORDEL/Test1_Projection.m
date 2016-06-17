iterMax = 2;

resGibbs = nan(3, iterMax);
iterSum = nan(3, iterMax);
distances = nan(3, iterMax);
feasibilityError = nan(3, iterMax);

% random phi
for i=1:iterMax
    phiInfeasible = rand(size(phi));
    phiStart = rand(size(phi));
    lambdaStart = rand(size(phi));
    
    [phi2,iterSum(1,i),~,resGibbs(1,i)] = projection2Gibbs(phiInfeasible,phiStart,matrices,lambdaStart,TriInfo);
    feasibilityError(1,i) = max(max(abs(sum(phi2,2)-1)), max(-min(phi2(:)),0));    
end
% feasible phi without zeros
for i=1:iterMax
    phiFeasiblePart = 0.2*rand(size(phi,1), TriInfo.sizePhi-1);
    phiFeasible = [phiFeasiblePart 1-sum(phiFeasiblePart,2)];
    phiStart = rand(size(phi));
    lambdaStart = rand(size(phi));
    [phi2,iterSum(2,i),~,resGibbs(2,i)] = projection2Gibbs(phiFeasible,phiStart,matrices,lambdaStart,TriInfo);
    feasibilityError(2,i) = max(max(abs(sum(phi2,2)-1)), max(-min(phi2(:)),0));
    
    phiDiff2 = phi2 - phiFeasible;

    
    
    distances(2,i) = norm(phiDiff2(:));

end
% feasible phi with zeros
for i=1:iterMax
    phiFeasiblePart = 0.2*rand(size(phi,1), TriInfo.sizePhi-1);
    phiFeasiblePartZeros = randi(4, size(phi,1), TriInfo.sizePhi-1) <= 2;
    phiFeasiblePart(phiFeasiblePartZeros) = 0;
    phiFeasible = [phiFeasiblePart 1-sum(phiFeasiblePart,2)];
    phiStart = rand(size(phi));
    lambdaStart = rand(size(phi));
    [phi2,iterSum(3,i),~,resGibbs(3,i)] = projection2Gibbs(phiFeasible,phiStart,matrices,lambdaStart,TriInfo);
    feasibilityError(3,i) = max(max(abs(sum(phi2,2)-1)), max(-min(phi2(:)),0));
    

    phiDiff2 = phi2 - phiFeasible;
%     distances(3,i,1) = sqrt(phiDiff1(:)'*matrices.H1scal*phiDiff1(:));

    
    distances(3,i) = norm(phiDiff2(:));
    
    
end
resGibbs
iterSum
feasibilityError
distances













