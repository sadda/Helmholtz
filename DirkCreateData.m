clear all;
close all;

addpath(genpath('./P1AFEM'));

%% Compute strains for the obtained solution

dirBase = fullfile('Results', 'Res_M1_0.000036_0.200', 'Results_Ref5_18_Elas');

load(fullfile(dirBase, 'DataAll.mat'));

[~, ~, ~, ~, ~,~, ~, ~, u, Theta, ~] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, dataEigen);
u   = u - material.eps0*[TriInfo.x; TriInfo.y];

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

% uxy(uxy>0) = 1;
% uxy(uxy<0) = -1;
%
% uyx(uyx>0) = 1;
% uyx(uyx<0) = -1;

PlotFunction(ux, TriInfo);
colormap gray;
PlotFunction(uxx, TriInfo);
colormap gray;
PlotFunction(uxy, TriInfo);
colormap gray;
PlotFunction(uy, TriInfo);
colormap gray;
PlotFunction(uyx, TriInfo);
colormap gray;
PlotFunction(uyy, TriInfo);
colormap gray;

phi         = ProlongPhi(phi, TriInfo);
coordinates = [TriInfo.x TriInfo.y];
elements    = TriInfo.e2p;

save('OT1_Data_Ela+Helm1', 'coordinates', 'elements', 'material', 'M', 'phi', 'ux', 'uy', 'uxx', 'uxy', 'uyx', 'uyy');

%% Generate the new mesh as Dirk wanted

requiredX      = TriInfo.requiredX;
requiredY      = TriInfo.requiredY;
requiredX(:,1) = [-0.75; 0.75];
data           = load('MeshesCreated/Mesh_New2/Data.mat');
TriInfoSmall   = data.TriInfo;

meshMaxElement = [max(diff(unique(TriInfoSmall.x))); max(diff(unique(TriInfoSmall.y)))];
epsilon        = 2*max(meshMaxElement);

mesh         = struct('coordinates', [TriInfoSmall.x TriInfoSmall.y], 'elements', TriInfoSmall.e2p);
coordinates  = mesh.coordinates;
elements     = mesh.elements;
dirichlet    = getBoundary(mesh);
neumann      = [];

refinementTol = 1e-6;
for i=1:5
    epsilon        = epsilon / 2;
    phi            = zeros(size(coordinates,1), 4);
    phi(:,end)     = 1;
    for j=1:size(phi,2)
        prescribed         = coordinates(:,1) >= requiredX(1,j) & coordinates(:,1) <= requiredX(2,j) & coordinates(:,2) >= requiredY(1,j) & coordinates(:,2) <= requiredY(2,j);
        phi(prescribed,:)  = 0;
        phi(prescribed,j)  = 1;
    end
    
    markedTriangles = zeros(size(elements,1), 1);
    for k = 1:size(elements,1)
        values = phi(elements(k,:),:);
        if any(values(:) < 1 - refinementTol & values(:) > refinementTol) || norm(values(1,:) - values(2,:)) > refinementTol || norm(values(1,:) - values(3,:)) > refinementTol
            markedTriangles(k) = 1;
        end
    end
    [coordinates, elements, ~, dirichlet, neumann] = refineNVBModified(coordinates, elements, [], dirichlet, neumann, logical(markedTriangles));
end

[TriInfo, Transformation]           = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),elements,requiredX,requiredY,TriInfo.sizePhi);
[Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
matrices.H1scalProlong              = TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;

width     = 0.01;
interface = zeros(size(coordinates,1), 1);
for j=1:4
    for i=1:2
        con1      = coordinates(:,1) >= requiredX(i,j)-width & coordinates(:,1) <= requiredX(i,j)+width & coordinates(:,2) >= requiredY(1,j)-width & coordinates(:,2) <= requiredY(2,j)+width;
        con2      = coordinates(:,2) >= requiredY(i,j)-width & coordinates(:,2) <= requiredY(i,j)+width & coordinates(:,1) >= requiredX(1,j)-width & coordinates(:,1) <= requiredX(2,j)+width;
        interface = interface | con1 | con2;
    end
end
con3      = coordinates(:,1) > min(coordinates(:,1)) & coordinates(:,1) < max(coordinates(:,1)) & coordinates(:,2) > min(coordinates(:,2)) & coordinates(:,2) < max(coordinates(:,2));
con4      = coordinates(:,2) ~= requiredY(1,4);
interface = interface & con3 & con4;

phi            = zeros(size(coordinates,1), 4);
phi(:,end)     = 1;
for j=1:size(phi,2)
    prescribed         = coordinates(:,1) >= requiredX(1,j) & coordinates(:,1) <= requiredX(2,j) & coordinates(:,2) >= requiredY(1,j) & coordinates(:,2) <= requiredY(2,j);
    phi(prescribed,:)  = 0;
    phi(prescribed,j)  = 1;
end
phiNew = zeros(size(coordinates,1), 4);
for j=1:4
    phiNew(:,j) = griddata(coordinates(~interface,1), coordinates(~interface,2), phi(~interface,j), coordinates(:,1), coordinates(:,2));
end

wavelength            = 1.64;
[constants, material] = ObtainData(epsilon, constants.alpha, wavelength, options);

PlotFunction(phiNew, TriInfo);

%% Compute the same as in the begining

options.method              = 0;
phi(~TriInfo.phiRowsFree,:) = [];
[~, ~, ~, ~, ~,~, ~, ~, u, Theta, ~] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, dataEigen);
u   = u - material.eps0*[TriInfo.x; TriInfo.y];

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

% uxy(uxy>0) = 1;
% uxy(uxy<0) = -1;
%
% uyx(uyx>0) = 1;
% uyx(uyx<0) = -1;

PlotFunction(ux, TriInfo);
colormap gray;
PlotFunction(uxx, TriInfo);
colormap gray;
PlotFunction(uxy, TriInfo);
colormap gray;
PlotFunction(uy, TriInfo);
colormap gray;
PlotFunction(uyx, TriInfo);
colormap gray;
PlotFunction(uyy, TriInfo);
colormap gray;

phi         = ProlongPhi(phi, TriInfo);
coordinates = [TriInfo.x TriInfo.y];
elements    = TriInfo.e2p;

save('OT1_Data_Ela+Helm2', 'coordinates', 'elements', 'material', 'M', 'phi', 'ux', 'uy', 'uxx', 'uxy', 'uyx', 'uyy');

%% Add show the prescribed values

fig = PlotFunction(phi, TriInfo);
hold on;
plot3([-2 2], [2.5 2.5], [4 4], 'k');

saveas(fig, 'Prescribed.jpg');
CropImage('Prescribed.jpg', 'Prescribed.jpg');










