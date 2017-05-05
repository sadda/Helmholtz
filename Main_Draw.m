close all;
clear all;

%% Prepare data

%folderNameBase = fullfile('Results', 'Res_M0_0.000400_2_1');
folderNameBase = fullfile('Res_M0_0.000400_2_1');
refineMesh     = 6;
drawPrescribed = 0;

%% Iteration count

iterCount  = zeros(1, refineMesh);
nodesCount = zeros(1, refineMesh);
resCount   = zeros(1, refineMesh);
glCount    = zeros(1, refineMesh);
for j = 1:refineMesh
    folderName = fullfile(folderNameBase, sprintf('Results_Ref%d_%d_Elas', j-1, 0));
    fileName   = fullfile(folderName, 'DataAll.mat');
    load(fileName, 'resAll', 'phi', 'TriInfo', 'matrices', 'constants');
    
    iter  = length(resAll(~isnan(resAll)));
    phiPr = ProlongPhi(phi, TriInfo);
    J2    = 0.5*(matrices.GradSq*phiPr(:))'*phiPr(:);
    J3    = 0.5*(matrices.Id'*phiPr(:) - (matrices.Mloc*phiPr(:))'*phiPr(:));
    
    iterCount(j)  = iter;
    nodesCount(j) = size(phi,1);
    resCount(j)   = min(resAll);
    glCount(j)    = constants.epsilon*J2 + 1/constants.epsilon*J3;
end

%% Draw phi and Theta

wavelength2 = 4;

folderName = fullfile(folderNameBase, sprintf('Results_Ref%d_%d_Elas', refineMesh-1, 0));
load(fullfile(folderName, 'DataAll.mat'), 'phi', 'TriInfo', 'Transformation', 'matrices', 'constants', 'material', 'options');
options.method = 0;
[~, ~, ~, ~, ~, ~, ~, ~, ~, Theta1] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
[~, material] = ObtainData(constants.epsilon, constants.alpha, wavelength2, options);
[~, ~, ~, ~, ~, ~, ~, ~, ~, Theta2] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);

fileName = 'Phi2.jpg';
fig      = PlotFunction(phi, TriInfo);
saveas(fig, fileName);
CropImage(fileName, fileName);

fileName = 'Theta2_1.jpg';
fig      = PlotFunction(Theta1, TriInfo);
saveas(fig, fileName);
CropImage(fileName, fileName);

fileName = 'Theta2_2.jpg';
fig      = PlotFunction(Theta2, TriInfo);
saveas(fig, fileName);
CropImage(fileName, fileName);

%% Plot mesh

for j = 1:refineMesh
    folderName = fullfile(folderNameBase, sprintf('Results_Ref%d_0_Elas', j-1));
    load(fullfile(folderName, 'DataAll.mat'), 'TriInfo');
    fig = figure;
    patch(TriInfo.x(TriInfo.e2p)', TriInfo.y(TriInfo.e2p)', 'w');
    set(gca,'XTickLabel','','YTickLabel','');
    
    fileName = sprintf('Mesh_%d.jpg', j-1);
    saveas(fig, fileName);
    CropImage(fileName, fileName);
end

%% Draw prescribed values

if drawPrescribed
    addpath(genpath('./P1AFEM'));
    
    data            = load('MeshesCreated/Mesh_New2/Data.mat');
    TriInfo0        = data.TriInfo;
    
    mesh          = struct('coordinates', [TriInfo0.x TriInfo0.y], 'elements', TriInfo0.e2p);
    coordinates   = mesh.coordinates;
    elements      = mesh.elements;
    dirichlet     = getBoundary(mesh);
    neumann       = [];
    refinementTol = 1e-6;
    for i=1:7
        markedTriangles = zeros(size(elements,1), 1);
        
        phi0        = zeros(size(coordinates,1), TriInfo0.sizePhi);
        phi0(:,end) = 1;
        for j=1:size(phi0, 1)
            for k=1:TriInfo0.sizePhi-1
                if coordinates(j,1) >= TriInfo0.requiredX(1,k) && coordinates(j,1) <= TriInfo0.requiredX(2,k) && coordinates(j,2) >= TriInfo0.requiredY(1,k) &&coordinates(j,2) <= TriInfo0.requiredY(2,k)
                    phi0(j,k)   = 1;
                    phi0(j,end) = 0;
                end
            end
        end
        for k = 1:size(elements,1)
            values = phi0(elements(k,:),:);
            % Mark all triangles in the interfacial interface
            if any(values(:) < 1 - refinementTol & values(:) > refinementTol) || norm(values(1,:) - values(2,:)) > refinementTol || norm(values(1,:) - values(3,:)) > refinementTol
                side1 = max(diff(sort(coordinates(elements(k,:),1))));
                side2 = max(diff(sort(coordinates(elements(k,:),2))));
                markedTriangles(k) = 1;
            end
        end
        [coordinates, elements, phi0, dirichlet, neumann] = refineNVBModified(coordinates, elements, phi0, dirichlet, neumann, logical(markedTriangles));
    end
    [TriInfo0, Transformation0]            = MeshCreateMatrices1(coordinates(:,1),coordinates(:,2),elements,TriInfo0.requiredX,TriInfo0.requiredY,TriInfo0.sizePhi);
    
    fileName = 'Prescribed.jpg';
    fig = PlotFunction(phi0, TriInfo0);
    hold on;
    plot3([min(TriInfo.x) max(TriInfo.x)], [TriInfo0.requiredY(1,end) TriInfo0.requiredY(1,end)], [10 10], 'k');
    saveas(fig, fileName);
    CropImage(fileName, fileName);
end

data      = [iterCount; nodesCount; resCount; glCount];
formatAll = {'%d', '%d', '%1.2e', '%1.3f'};
TableToTex(data, formatAll);

%% Compute perimeter

perAll = zeros(1, refineMesh);
for j = 1:refineMesh
    folderName = fullfile(folderNameBase, sprintf('Results_Ref%d_%d_Elas', j-1, 0));
    load(fullfile(folderName, 'DataAll.mat'), 'phi', 'TriInfo');
    
    x        = TriInfo.x;
    y        = TriInfo.y;
    e2p      = TriInfo.e2p;
    nelement = TriInfo.nelement;
    
    phiPr    = ProlongPhi(phi, TriInfo);
    phiMax   = zeros(size(phiPr));
    [~,phiI] = max(phiPr, [], 2);
    for i=1:size(phiMax,1)
        phiMax(i,phiI(i)) = 1;
    end
    
    mesh = struct('coordinates', [TriInfo.x TriInfo.y], 'elements', TriInfo.e2p);
    element2neighbors = getNeighbors(mesh);
    
    per = 0;
    figure;
    hold on;
    
    for k=1:nelement
        neighbors = element2neighbors(k,:);
        for i=1:3
            neighbor = neighbors(i);
            if neighbor > 0
                
                value1 = phiMax(e2p(k,:),:);
                value2 = phiMax(e2p(neighbor,:),:);
                
                if ~isequal(value1(1,:), value1(2,:)) || ~isequal(value1(1,:), value1(3,:)) || ~isequal(value2(1,:), value2(2,:)) || ~isequal(value2(1,:), value2(3,:))
                    
                    int    = intersect(e2p(k,:), e2p(neighbor,:));
                    value3 = phiMax(int,:);
                    if isequal(value3(1,:), value3(2,:))
                        x1  = x(int);
                        y1  = y(int);
                        per = per + norm([x1(1)-x1(2), y1(1)-y1(2)]);
                        
                        line(x1, y1);
                    end
                end
            end
        end
    end
    perAll(j) = per / 4;
end

TableToTex(pi/4*perAll, {'%1.3f'});

