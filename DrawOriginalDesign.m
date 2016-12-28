addpath('OldCodes');

load('MeshesCreated/Mesh_Big/Data.mat');

phi1  = TriInfo.x >= -0.75 & TriInfo.x <= 0.75 & TriInfo.y >= 1.1 & TriInfo.y <= 1.4;
phi2  = TriInfo.x >= -0.75 & TriInfo.x <= 0.75 & TriInfo.y > 1.4 & TriInfo.y <= 1.8;
phi3  = TriInfo.y < 1.1;
phi4  = ~(phi1 | phi2 | phi3);
fixGe = phi1;

phi         = zeros(TriInfo.npoint, TriInfo.sizePhi);
phi(phi1,1) = 1;
phi(phi2,2) = 1;
phi(phi3,3) = 1;
phi(phi4,4) = 1;

[TriInfo, matrices, phi]  = ModifyMatrices(fixGe, TriInfo, matrices, Transformation, phi);
[constants, material]     = ObtainData(1e-5, 0, 2*pi);
u                         = ComputeElasticity(phi,TriInfo,Transformation,matrices,constants,material);
[J, data.J1, J2, J3]      = ComputeObjective(phi,u,TriInfo,matrices,constants);


options                   = struct('computeG', 0, 'computeU', 1, 'symmetrize', 0, 'separateObjective', 0);
[~,~,~,~,~,~,~,~,~,Theta] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options);


totalArea    = TriInfo.totalArea;
phiProlonged = ProlongPhi(phi, TriInfo);


minX      = min(TriInfo.x);
maxX      = max(TriInfo.x);
minY      = min(TriInfo.y);
maxY      = max(TriInfo.y);


% fig = figure();
% hold on;
% trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, phiProlonged*(1:TriInfo.sizePhi)');
% plot([minX maxX maxX minX minX], [minY minY maxY maxY minY], 'k');
% title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', totalArea, data.J1, data.J1 / totalArea));
% view(2);
% shading interp;
% colormap(colormap(gray));
% set(gca,'xcolor',get(gcf,'color'));
% set(gca,'xtick',[]);
% set(gca,'ycolor',get(gcf,'color'));
% set(gca,'ytick',[]);
% saveas(fig, 'Or1a.jpg');






fig = figure();
trisurf(TriInfo.e2p, TriInfo.x, TriInfo.y, Theta);
title(sprintf('Area = %1.5f, strain = %1.5f, ratio = %1.5f', totalArea, data.J1, data.J1 / totalArea));
view(2);
shading interp;
colormap jet;
colorbar;
hold on;
set(gca,'xcolor',get(gcf,'color'));
set(gca,'xtick',[]);
set(gca,'ycolor',get(gcf,'color'));
set(gca,'ytick',[]);
% Draw boundary edges for fixGePrev
fixGe    = TriInfo.phiProlongationVector(:,1) == 1;
xEle     = TriInfo.x(TriInfo.e2p);
yEle     = TriInfo.y(TriInfo.e2p);
fixGeEle = sum(fixGe(TriInfo.e2p), 2) == 3;
xDraw = xEle(fixGeEle,:);
xDraw = xDraw(:,[1 1 2 2 3 3]);
xDraw = reshape(xDraw, [], 2);
yDraw = yEle(fixGeEle,:);
yDraw = yDraw(:,[1 1 2 2 3 3]);
yDraw = reshape(yDraw, [], 2);
sortIndex          = (xDraw(:,1) > xDraw(:,2)) | (xDraw(:,1) == xDraw(:,2) & yDraw(:,1) > yDraw(:,2));
xDraw(sortIndex,:) = xDraw(sortIndex,[2 1]);
yDraw(sortIndex,:) = yDraw(sortIndex,[2 1]);
xyDraw   = sortrows([xDraw yDraw]);
xyDiff1  = diff(xyDraw);
xyDiff2  = [ones(1, 4); xyDiff1];
xyDiff3  = [xyDiff1; ones(1, 4)];
xyUnique = sum(abs(xyDiff2),2)~=0 & sum(abs(xyDiff3),2)~=0;
if ~isempty(xyDraw)
    xyDraw   = xyDraw(xyUnique,:);
    plot3(xyDraw(:,1:2)', xyDraw(:,3:4)', 2*max(Theta)*ones(size(xyDraw,1),2)', 'k');
else
    scatter3(TriInfo.x(fixGe), TriInfo.y(fixGe), 2*max(Theta)*ones(sum(fixGe),1), 'k');
end
saveas(fig, 'Or_6.28.jpg')




