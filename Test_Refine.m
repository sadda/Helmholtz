clear all;
close all;

%% Load and prepare data

addpath(genpath('./P1AFEM'));
load('MeshesCreated/Mesh_GeFree/Data.mat');

mesh.elements    = TriInfo.e2p;
mesh.coordinates = [TriInfo.x TriInfo.y];
dirichlet        = getBoundary(mesh);
neumann          = [];
npoint0          = TriInfo.npoint;
mesh             = genBisectionMesh(mesh);

figure;
for k=1:size(mesh.elements,1)
    patch(mesh.coordinates(mesh.elements(k,:),1)',mesh.coordinates(mesh.elements(k,:),2)',1);
    hold on;
end

%% Refine the mesh

markedTriangles = true(size(mesh.elements,1), 1);
[coordinates, elements, ~, dirichlet, neumann] = refineNVBModified([TriInfo.x TriInfo.y], TriInfo.e2p, [], dirichlet, neumann, logical(markedTriangles));

figure;
for k=1:size(elements,1)
    patch(coordinates(elements(k,:),1)',coordinates(elements(k,:),2)',1);
    hold on;
end

%% Coarsen the middle part of the mesh

markedTriangles = false(size(elements,1), 1);
y = coordinates(:,2);
marked = sum(y(elements) >= 1 & y(elements) <= 2, 2) == 3;
markedTriangles(marked) = true;
[coordinates, elements, ~, dirichlet, neumann] = coarsenNVBModified(npoint0, coordinates, elements, [], dirichlet, neumann, logical(markedTriangles));

figure;
for k=1:size(elements,1)
    patch(coordinates(elements(k,:),1)',coordinates(elements(k,:),2)',1);
    hold on;
end



