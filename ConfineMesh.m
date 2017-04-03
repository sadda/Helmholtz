clear all;
close all;

addpath(genpath('./P1AFEM'));

load(fullfile('Res', 'DataAll2.mat'), 'TriInfo', 'matrices', 'phi');

%% Gather data

phi           = ProlongPhi(phi, TriInfo);
TriInfoOld    = TriInfo;
phiOld        = phi;
mesh          = struct('coordinates', [TriInfoOld.x TriInfoOld.y], 'elements', TriInfoOld.e2p);
coordinates   = mesh.coordinates;
elements      = mesh.elements;
mesh          = struct('coordinates', coordinates, 'elements', elements);
dirichlet     = getBoundary(mesh);
neumann       = [];
dataToProlong = phi;

%% Refine

markedTriangles = zeros(size(elements,1), 1);
for k = 1:size(elements,1)
    values = phi(elements(k,:),:);
    % Add some condition for element refinement
    if all(coordinates(elements(k,:),1) <= 0)
        markedTriangles(k) = 1;
    end
end
[coordinates, elements, dataToProlong, dirichlet, neumann] = refineNVBModified(coordinates, elements, dataToProlong, dirichlet, neumann, logical(markedTriangles));
phi = dataToProlong;

%% Plot old mesh

figure;
patch(TriInfoOld.x(TriInfoOld.e2p)', TriInfoOld.y(TriInfoOld.e2p)', 'w');

figure;
trisurf(TriInfoOld.e2p, TriInfoOld.x, TriInfo.y, phiOld*(1:4)');
view(2);
shading interp;

%% Plot new mesh

x = coordinates(:,1);
y = coordinates(:,2);
figure;
patch(x(elements)', y(elements)', 'w');

figure;
trisurf(elements, x, y, phi*(1:4)');
view(2);
shading interp;





