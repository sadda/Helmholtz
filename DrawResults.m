clear all;
close all;

J   = zeros(4,1);
J1  = zeros(4,1);
J2  = zeros(4,1);
J3  = zeros(4,1);
DoF = zeros(4,1);

for i=1:4
    switch i
        case 1
            load('Results/OT1_Elasticity_Problem2_1/Results_Ref4/DataAll.mat', 'TriInfo', 'Transformation', 'matrices', 'constants', 'material', 'phi');
            epsilon = constants.epsilon;
            alpha   = constants.alpha;
        case 2
            load('Results/OT1_Elasticity_Problem2_2/InteriorPoint/Sol5.mat', 'TriInfo', 'Transformation', 'matrices', 'constants', 'material', 'phi');
        case 3
            load('Results/OT1_Elasticity_Problem2_3/InteriorPoint2/Sol5.mat', 'TriInfo', 'Transformation', 'matrices', 'constants', 'material', 'phi');
        case 4
            load('MeshesCreated/Void2/DataAll.mat');
            load('MeshesCreated/Void2/phi.mat');
            load('Results/OT1_Elasticity_Problem2_3/InteriorPoint2/Sol5.mat', 'constants');
            epsilon = constants.epsilon;
            alpha   = constants.alpha;
            phi(:,4)     = phi(:,4) + phi(:,2);
            phi(:,5)     = phi(:,5) + phi(:,3);
            phi(:,[2 3]) = [];
            TriInfo.sizePhi = 4;
            TriInfo.phiProlongationVector = zeros(TriInfo.npoint, TriInfo.sizePhi);
            [constants, material] = ObtainData(epsilon, alpha);
            u     = elasticity_adjoint(phi,TriInfo,Transformation, matrices,constants,material);
            J1(i) = -(matrices.TrXD+matrices.TrYD)'*u;
    end
    %% Compute objective
    DoF(i) = size(phi,1);
    if i < 4
        u = elasticity_adjoint(phi,TriInfo,Transformation, matrices,constants,material);
        [J(i), J1(i), J2(i), J3(i)] = ComputeObjective(phi,u,TriInfo,matrices,constants);
    end
    
    %% Draw phases
    phi  = ProlongPhi(phi, TriInfo);
    
    minX = min(TriInfo.x);
    maxX = max(TriInfo.x);
    minY = min(TriInfo.y);
    maxY = max(TriInfo.y);
    
    fig = figure;
    hold on;
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phi*(1:TriInfo.sizePhi)');
    plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
    view(2);
    shading interp;
    colormap(colormap(gray));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    saveas(fig, ['Microbridge', int2str(i), '.jpg']);
    
    %% Draw biaxial strain
    v      = matrices.Mloc2D\(matrices.Tr2D*u);
    vx     = v(1:TriInfo.npoint);
    vy     = v(TriInfo.npoint+1:2*TriInfo.npoint);
    tr_eps = abs(vx + vy)/2;
    fig = figure;
    hold on;
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, tr_eps);
    plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
    view(2);
    shading interp;
    colormap(flipud(colormap(gray)));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    colorbar;
    caxis([0 1.4e-2]);
    saveas(fig, ['Microbridge_BiaxialStrain', int2str(i), '.jpg']);
end

(J1./J1(4)-1)*100


%%

nodesNumber = zeros(3,5);
for i=1:3
    for j=1:5
        switch i
            case 1
                load(['Results/OT1_Elasticity_Problem2_1/Results_Ref', int2str(j-1), '/DataAll.mat'], 'TriInfo');
            case 2
                load(['Results/OT1_Elasticity_Problem2_2/InteriorPoint/Sol', int2str(j), '.mat'], 'TriInfo');
            case 3
                load(['Results/OT1_Elasticity_Problem2_3/InteriorPoint2/Sol', int2str(j), '.mat'], 'TriInfo');
        end
        nodesNumber(i,j) = TriInfo.npoint;
    end
end


%%


for i=1:3
    switch i
        case 1
            data   = dlmread('Results/OT1_Elasticity_Problem2_1/Results_Ref3/Iteration_data.csv', ' ', 1, 1);
            resAll = data(:,4);
        case 2
            load('Results/OT1_Elasticity_Problem2_2/InteriorPoint/Sol4.mat', 'resAll');
        case 3
            load('Results/OT1_Elasticity_Problem2_3/InteriorPoint2/Sol4.mat', 'resAll');
    end
    resAll = resAll(resAll > 0);
    min(resAll)
    fig = figure;
    plot(log10(resAll), 'k');
    xlim([1 120]);
    ylim([-10 1]);
    xlabel('iteration');
    ylabel('log_{10}(residual)');
    saveas(fig, ['Residual', int2str(i), '.jpg']);
end



%%


for i=1:2
    switch i
        case 1
            load('Results/OT1_Elasticity_Problem2_5/Results_Ref3/DataAll.mat', 'TriInfo', 'phi');
        case 2
            load('Results/OT1_Elasticity_Problem2_6/Results_Ref3/DataAll.mat', 'TriInfo', 'phi');
    end
    
    phi  = ProlongPhi(phi, TriInfo);
    
    minX = min(TriInfo.x);
    maxX = max(TriInfo.x);
    minY = min(TriInfo.y);
    maxY = max(TriInfo.y);
    
    fig = figure;
    hold on;
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phi*(1:TriInfo.sizePhi)');
    plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
    view(2);
    shading interp;
    colormap(colormap(gray));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    saveas(fig, ['MicrobridgeScaled', int2str(i), '.jpg']);
end


%%

alphaAll = 10.^(linspace(log10(1e-6),log10(1e-2),48));
load(fullfile('MeshesCreated', 'VoidUR1Pr', 'Data.mat'), 'TriInfo');
for i=1:length(alphaAll)
    alpha   = alphaAll(i);
    load(fullfile('Results', 'OT1_Elasticity_Problem2_4', 'RES', sprintf('Data_Alpha=%1.7f.mat', alpha)), 'phi');
    
    phi  = ProlongPhi(phi, TriInfo);
    
    minX = min(TriInfo.x);
    maxX = max(TriInfo.x);
    minY = min(TriInfo.y);
    maxY = max(TriInfo.y);
    
    fig = figure('visible', 'off');
    hold on;
    trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phi*(1:TriInfo.sizePhi)');
    plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
    view(2);
    shading interp;
    colormap(colormap(gray));
    set(gca,'xcolor',get(gcf,'color'));
    set(gca,'xtick',[]);
    set(gca,'ycolor',get(gcf,'color'));
    set(gca,'ytick',[]);
    saveas(fig, fullfile('Results', 'OT1_Elasticity_Problem2_4', 'RES', sprintf('Shape2_Alpha=%1.7f.jpg', alphaAll(i))));
end

%%

alphaAll = 10.^(linspace(log10(1e-6),log10(1e-2),48));
J1All    = zeros(size(alphaAll));
for i=1:length(alphaAll);
    alpha   = alphaAll(i);
    load(fullfile('Results', 'OT1_Elasticity_Problem2_4', 'RES', sprintf('Data_Alpha=%1.7f.mat', alpha)), 'J1');
    J1All(i) = J1;
end
fig = figure;
% plot(log10(alphaAll), -J1All, '-*k');
plot(log10(alphaAll), -J1All, '-k');
hold on;
%line([log10(2e-4) log10(2e-4)], [min(-J1All) max(-J1All)]);
xlabel('log_{10}(\alpha)');
ylabel('strain');
saveas(fig, 'StrainDependence.jpg');



