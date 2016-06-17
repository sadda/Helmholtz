clear all;
close all;

meshName = 'MeshesCreated/VoidUR1SymModPr/';

% Start with the original solution
load([meshName, '/Data.mat']);
meshMaxElement = [max(diff(unique(TriInfo.x))); max(diff(unique(TriInfo.y)))];
epsilon = 2*max(meshMaxElement);
load('MeshesCreated/Void2/phiInit.mat');
load('MeshesCreated/Void2/TriInfo.mat');
TriInfoOriginal = TriInfo;
load([meshName, '/Data.mat']);
phi = ProlongPhiMesh(phi, TriInfoOriginal, TriInfo);
phi(~TriInfo.phiRowsFree,:) = [];

% Prolong the solution
TriInfoOld = TriInfo;
phiProlonged = ProlongPhi(phi, TriInfoOld);
[Transformation, TriInfo, matrices] = MeshCreateProlong(TriInfoOld, phiProlonged, 0, 1.5*meshMaxElement);
meshMaxElement = meshMaxElement / 2;
phi = ProlongPhiMesh(phiProlonged, TriInfoOld, TriInfo);
phi(~TriInfo.phiRowsFree,:) = [];
epsilon = epsilon / 2;







% Assign the parameters
[constants, material] = ObtainData(epsilon, 1, 1, TriInfo.sizePhi);
dirNameBase           = ['Results_Ref', int2str(meshIndex-1), '_ParUpdate_'];
dirNameIndex          = 0;
try
    save(['Ref', int2str(meshIndex)]);
end
% Run the optimization
if meshIndex == 1
    for i=1:refineAlpha
        [u,p]               = elasticity_adjoint(phi,TriInfo,Transformation,matrices,constants,material);
        [J, J1, J2, J3, J4] = functionValue(phi,u,TriInfo,matrices,constants);
        alpha               = abs(J1) / (epsilon*J2 + 1/epsilon*J3);
        alpha               = 1/2*alpha;
        alpha               = 5^(refineAlpha-i)*alpha;
        gamma               = 0;
        dirName             = [dirNameBase, int2str(dirNameIndex)];
        dirNameIndex        = dirNameIndex + 1;
        [phi, tInit]        = RunOptimization(alpha, epsilon, gamma, TriInfo, Transformation, matrices, dirName, IterMax, drawResults, phi, algorithm, tInit);
    end
else
    [u,p]               = elasticity_adjoint(phi,TriInfo,Transformation,matrices,constants,material);
    [J, J1, J2, J3, J4] = functionValue(phi,u,TriInfo,matrices,constants);
    alpha               = abs(J1) / (epsilon*J2 + 1/epsilon*J3);
    alpha               = 1/2*alpha;
    gamma               = 0;
    dirName             = [dirNameBase, int2str(dirNameIndex)];
    dirNameIndex        = dirNameIndex + 1;
    [phi, tInit]        = RunOptimization(alpha, epsilon, gamma, TriInfo, Transformation, matrices, dirName, IterMax, drawResults, phi, algorithm, tInit);
    if meshIndex == refineMesh
        [u,p]                  = elasticity_adjoint(phi,TriInfo,Transformation,matrices,constants,material);
        [~,J1,~,~,~,J4_1,J4_2] = functionValue(phi,u,TriInfo,matrices,constants);
        gamma0                 = abs(J1) / (J4_1 + J4_2);
        gammaAll               = 1/100*gamma0*gammaBase.^(0:gammaCount-1);
        for gamma=gammaAll
            dirName            = [dirNameBase, int2str(dirNameIndex)];
            dirNameIndex       = dirNameIndex + 1;
            [phi, tInit]       = RunOptimization(alpha, epsilon, [gamma; gamma], TriInfo, Transformation, matrices, dirName, IterMax, drawResults, phi, algorithm, tInit);
        end
    end
end
end

% exit;



