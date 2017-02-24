clear all;
close all;

addpath('./OldCodes');
addpath(genpath('./P1AFEM'));

refineMesh    = 4;
drawResults   = 0;
iterMax       = 25;
iterMaxIn     = 50;

pars                       = [];
%         pars.alpha                 = [1e-5 5e-5 1e-4 5e-4 1e-3];
%         pars.jointObjectiveThetaLp = 2;
%         pars.jointObjectivePhiLp   = 2;
%         pars.geLevel2              = 1;
%         pars.gammaPen              = [1e-2 5e-2 1e-1 5e-1 1e-0];
%         pars.alphaGe               = [1e-3 5e-3 1e-2 5e-2 1e-1];
pars.alpha                 = 1e-4;
pars.jointObjectiveThetaLp = 2;
pars.jointObjectivePhiLp   = 2;
pars.geLevel2              = 1;
pars.gammaPen              = 0.1;
pars.alphaGe               = 5e-2;
parsAll                    = combvec(unique(pars.alpha), pars.jointObjectiveThetaLp, pars.jointObjectivePhiLp, pars.geLevel2, unique(pars.gammaPen), unique(pars.alphaGe));


% parfor parsIndex = 1:size(parsAll, 2)
for parsIndex = 1:size(parsAll, 2)
    
    meshIndex = refineMesh;
    alpha     = parsAll(1,parsIndex);
    tInit1 = 1e3;
    tInit2 = 1e3;
    
    options.computeU              = 1;
    options.separateObjective     = 0;
    options.jointObjectiveThetaLp = parsAll(2,parsIndex);
    options.jointObjectivePhiLp   = parsAll(3,parsIndex);
    options.normalizationLp       = 2;
    options.geLevel2              = parsAll(4,parsIndex);
    options.gammaPen              = parsAll(5,parsIndex);
    options.alphaGe               = parsAll(6,parsIndex);
    
    dirNameBase = sprintf('Res_M0_%1.6f_%d_%d_%1.4f_%1.4f_%1.4f', alpha, options.jointObjectiveThetaLp, options.jointObjectivePhiLp, options.geLevel2, options.gammaPen, options.alphaGe);
    
    data           = load(fullfile('qwe', dirNameBase, 'Results_Ref3_0_Elas', 'DataAll.mat'));
    phi            = data.phi;
    TriInfo        = data.TriInfo;
    Transformation = data.Transformation;
    matrices       = data.matrices;
    material       = data.material;
    constants      = data.constants;
    options        = data.options;
    tInit1         = data.t;
    
    
    sum(phi)
    
    %     options.alphaGe  = 0;
    options.gammaPen = 10*options.gammaPen;
    
    
    iterIn  = 0;
    dirName     = sprintf('MOD_Results_Ref%d_%d_Elas', meshIndex-1, iterIn);
    
    picName        = fullfile(dirName, 'Pictures');
    if exist(picName, 'dir')
        rmdir(picName, 's');
    end
    mkdir(picName);
    
    
    [phi, tInit1]             = ProjectedGradients(TriInfo, Transformation, matrices, material, constants, dirName, iterMax, drawResults, phi, tInit1, options);
    [~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
    
    figure;
    phiPr = ProlongPhi(phi, TriInfo);
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiPr(:,1));
    view(3);
    
    figure;
    phiPr = ProlongPhi(phi, TriInfo);
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiPr(:,2));
    view(3);
    
    data.J1        = NaN;
    data.totalArea = NaN;
    data.iterIn    = iterIn;
    data.meshIndex = meshIndex;
    data.drawFixGe = 0;
    
    DrawResults(phi, Theta, TriInfo, picName, data);
    
end




% exit;









