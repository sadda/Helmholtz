close all;
clear all;

%% Prepare data

folderName1 = fullfile('Results', 'Res_M1_0.000036_0.200');
folderName2 = fullfile('Results', 'Res_M0_0.000114_2_1_0.0000_0.0000_0.0000');

refineMesh  = 6;
folderCount = zeros(refineMesh, 1);

for i = 1:refineMesh
    j = 1;
    while true
        folderName = fullfile(folderName1, sprintf('Results_Ref%d_%d_Elas', i-1, j));
        if ~exist(folderName, 'dir')
            break;
        end
        j = j+1;
    end
    folderCount(i) = j-1;
end


%% Iteration count

iterCount  = zeros(2, refineMesh);
nodesCount = zeros(2, refineMesh);
for i = 1:refineMesh
    for j = 1:folderCount(i)
        folderName = fullfile(folderName1, sprintf('Results_Ref%d_%d_Elas', i-1, j));
        fileName   = fullfile(folderName, 'DataAll.mat');
        load(fileName, 'resAll');
        
        iter = length(resAll(~isnan(resAll)));
        iterCount(1,i) = iterCount(1,i) + iter;
    end
    load(fileName, 'phi');
    nodesCount(1,i) = size(phi,1);
    
    folderName = fullfile(folderName2, sprintf('Results_Ref%d_%d_Elas', i-1, 0));
    fileName   = fullfile(folderName, 'DataAll.mat');
    load(fileName, 'resAll', 'phi');
    
    iter = length(resAll(~isnan(resAll)));
    iterCount(2,i)  = iter;
    nodesCount(2,i) = size(phi,1);
end

%% Draw phi and Theta

wavelength2 = 4;
folderName  = fullfile(folderName1, sprintf('Results_Ref%d_%d_Elas', refineMesh-1, folderCount(refineMesh)));
load(fullfile(folderName, 'DataAll.mat'), 'phi', 'TriInfo', 'Transformation', 'matrices', 'constants', 'material', 'options');
options.method = 0;
[~, ~, ~, ~, ~, ~, ~, ~, ~, Theta1] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);
[~, material] = ObtainData(constants.epsilon, constants.alpha, wavelength2, options);
[~, ~, ~, ~, ~, ~, ~, ~, ~, Theta2] = ComputeData(phi, TriInfo, Transformation, matrices, constants, material, options, []);

fileName = 'Phi1.jpg';
fig      = PlotFunction(phi, TriInfo);
saveas(fig, fileName);
CropImage(fileName, fileName);

fileName = 'Theta1_1.jpg';
fig      = PlotFunction(Theta1, TriInfo);
saveas(fig, fileName);
CropImage(fileName, fileName);

fileName = 'Theta1_2.jpg';
fig      = PlotFunction(Theta2, TriInfo);
saveas(fig, fileName);
CropImage(fileName, fileName);


folderName = fullfile(folderName2, sprintf('Results_Ref%d_%d_Elas', refineMesh-1, 0));
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

for i = 1:refineMesh
    folderName = fullfile(folderName1, sprintf('Results_Ref%d_1_Elas', i-1));
    % folderName = fullfile(folderName2, sprintf('Results_Ref%d_0_Elas', i-1));
    load(fullfile(folderName, 'DataAll.mat'), 'TriInfo');
    fig = figure;
    patch(TriInfo.x(TriInfo.e2p)', TriInfo.y(TriInfo.e2p)', 'w');
    set(gca,'XTickLabel','','YTickLabel','');
    
    fileName = sprintf('Mesh_%d.jpg', i-1);
    saveas(fig, fileName);
    CropImage(fileName, fileName);
end




