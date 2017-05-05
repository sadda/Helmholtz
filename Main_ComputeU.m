clear all;
close all;

symmetrize = 0;

dirBase = fullfile('Results', 'Res_M0_0.000400_2_1', 'Results_Ref5_0_Elas');
load(fullfile(dirBase, 'DataAll.mat'), 'TriInfo', 'Transformation', 'matrices', 'constants', 'material', 'options', 'dataEigen');
load(fullfile(dirBase, 'DataAll.mat'), 'phi');

%% Symmetrize

phiPr  = ProlongPhi(phi, TriInfo);
if symmetrize
    phiSym = phiPr;
    for i=1:size(phiSym,1)
        if TriInfo.x(i) < 0
            ind = find(TriInfo.x == -TriInfo.x(i) & TriInfo.y == TriInfo.y(i));
            if length(ind) == 1
                phiSym(i,:) = phiPr(ind,:);
            else
                res = (TriInfo.x + TriInfo.x(i)).^2 + (TriInfo.y - TriInfo.y(i)).^2;
                [~,ind] = min(res);
                phiSym(i,:) = phiPr(ind,:);
            end
        end
    end
    phiSym = phiSym(TriInfo.phiRowsFree,:);
else
    phiSym = phi;
end

%% Compute data

[~, ~, ~, ~, ~,~, ~, ~, u, Theta, ~] = ComputeData(phiSym, TriInfo, Transformation, matrices, constants, material, options, dataEigen);

npoint = TriInfo.npoint;
M      = matrices.Mloc(1:npoint,1:npoint);

C   = matrices.Tr2D;
Cx  = C(1:npoint,1:npoint);
Cy  = C(npoint+1:2*npoint,npoint+1:2*npoint);
ux  = u(1:npoint);
uy  = u(npoint+1:end);
uxx = M \ (Cx*ux);
uxy = M \ (Cy*ux);
uyx = M \ (Cx*uy);
uyy = M \ (Cy*uy);

phiPr       = ProlongPhi(phiSym, TriInfo);
coordinates = [TriInfo.x TriInfo.y];
elements    = TriInfo.e2p;


%% Compote max uxx

indexPhi        = find(phiPr(:,1) > 0.99);
[~,indexPhiMax] = max(uxx(indexPhi));
index           = indexPhi(indexPhiMax);

phi1 = phiPr(:,1);
red1 = sum(phi1(elements) >= 0.9, 2) == 3;
red2 = unique(elements(red1,:));

uxx(index)
phiPr(index,:)
[TriInfo.x(index), TriInfo.y(index)]

%% Plot data



if symmetrize
    PlotFunction(phi, TriInfo);
    fig = PlotFunction(phiSym, TriInfo);
    saveas(fig, 'Des1_Design.jpg');
else
    fig = PlotFunction(phi, TriInfo);
    saveas(fig, 'Des1_Design.jpg')    
end

fig = PlotFunction(uxx, TriInfo);
caxis([-0.02 0.02]);
saveas(fig, 'Des1_Uxx_Scale1.jpg');

fig = PlotFunction(uxx, TriInfo);
caxis([-0.01 0.01]);
saveas(fig, 'Des1_Uxx_Scale2.jpg');

fig = figure();
hold on;
trisurf(TriInfo.e2p(red1,:), TriInfo.x, TriInfo.y, uxx);
view(2);
shading interp;
set(gca,'xcolor',get(gcf,'color'));
set(gca,'xtick',[]);
set(gca,'ycolor',get(gcf,'color'));
set(gca,'ytick',[]);
colormap(flipud(colormap(hot)));
colorbar;
caxis([0 0.008]);
saveas(fig, 'Des1_Uxx_Germanium.jpg');

%% Save data

phi  = ProlongPhi(phi, TriInfo);
save('OT1_Data_Ela+Helm4', 'coordinates', 'elements', 'material', 'M', 'phi', 'ux', 'uy', 'uxx', 'uxy', 'uyx', 'uyy');





