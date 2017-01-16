pars                  = [];
pars.jointObjectiveLp = repmat([1 2 8 2 8], 1, 3);
pars.normalizationLp  = repmat([1 2 8 8 2], 1, 3);
pars.alpha            = kron(linspace(1e-5, 1e-4, 3), ones(1, 5));
iterIn = 0;

for parsIndex = 1:length(pars.alpha)
    dirNameBase = sprintf('Results_Tog_Obj=%d_Norm=%d_Alpha=%1.6f', pars.jointObjectiveLp(parsIndex), pars.normalizationLp(parsIndex), pars.alpha(parsIndex));
    dirNameBase = fullfile('RES', dirNameBase);
    
    for meshIndex=1:4
        dirName1                 = sprintf('Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
        dirName2                 = sprintf('Results_Ref%d_%d_Helm', meshIndex-1, iterIn);
        dirName1                 = fullfile(dirNameBase, dirName1);
        dirName2                 = fullfile(dirNameBase, dirName2);
        
        load(fullfile(dirName1, 'phi.mat'));
        load(fullfile(dirName1, 'DataAll.mat'), 'TriInfo', 'matrices');
        
        phiPr = ProlongPhi(phi, TriInfo);
        
        for i=1:4        
            sumPhi(meshIndex,i) = ones(1,TriInfo.npoint)*matrices.Mloc(1:TriInfo.npoint,1:TriInfo.npoint)*phiPr(:,i);
        end
        sumPhi(meshIndex,5) = sum(phi(:,1) > 0.5);
        sumPhi(meshIndex,6) = sum(phi(:,1) > 0.75);
        sumPhi(meshIndex,7) = sum(phi(:,1) > 0.99);
        sumPhi(meshIndex,8) = max(phi(:,1));
    end
    sumPhi
end