clear all;
close all;

tInit       = 1e3;
refineMesh  = 4;
drawResults = 1;
IterMax     = 10;
alpha       = 1e-5;

for meshIndex = 1:refineMesh
    % Construct mesh and determine the starting point
    if meshIndex == 1
%         load('MeshesCreated/VoidUR0/Data.mat');
        load('MeshesCreated/VoidUR0_SmallGe/Data.mat');
        meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))];
        epsilon        = 2*max(meshMaxElement);
        phi            = 1/TriInfo.sizePhi*ones(TriInfo.npointRed, TriInfo.sizePhi);
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
    
    [constants, material] = ObtainData(epsilon, alpha);
    [phi, tInit] = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, IterMax, drawResults, phi, tInit);
end

exit;



