clear all;
close all;

tInit       = 1;
refineMesh  = 4;
drawResults = 1;
IterMax     = 5000;
algorithm   = 1;
alpha       = 2e-4;

for meshIndex = 1:refineMesh
    % Construct mesh and determine the starting point
    if meshIndex == 1
        load('MeshesCreated/VoidUR0Pr/Data.mat');
        meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))];
        epsilon        = 2*max(meshMaxElement);
        phi            = rand(TriInfo.npointRed, TriInfo.sizePhi);
        phi            = phi ./ repmat(sum(phi,2), 1, TriInfo.sizePhi);
    else
        phiProlonged = ProlongPhi(phi, TriInfo);
        [Transformation, TriInfo, matrices, phi] = MeshCreateProlong(TriInfo, phiProlonged, 0, 1.5*meshMaxElement, phiProlonged);
        phi(~TriInfo.phiRowsFree,:) = [];
        epsilon = epsilon / 2;
        meshMaxElement = meshMaxElement / 2;
    end
    
    dirName = ['Results_Ref', int2str(meshIndex-1)];
    try
        save(['Ref', int2str(meshIndex)]);
    end
    [phi, tInit] = ProjectedGradients_RunOptimization(alpha, epsilon, TriInfo, Transformation, matrices, dirName, IterMax, drawResults, phi, tInit);
end

exit;



