clear all;
close all;

% Always check the mesh visually before working with it :)

plotResults      = 1;
prolongGermanium = 1;

xRange    = [-2 2];
yRange    = [0 3];
requiredX = [-0.75 -0.75 -2 -2; 0.75 0.75 2 2];
requiredY = [1.1 1.48 0 2.5; 1.4 1.8 0.96 3];

NX = 33;
NY = 25;
addedX = [];
addedY = [1.445];
results_file='MeshesCreated/VoidUR0';
sizePhi = 4;

if prolongGermanium
    opticalCavity = [requiredX(:,1) requiredY(:,1)];
    requiredX(1,1) = requiredX(1,1) - 0.01;
    requiredX(2,1) = requiredX(2,1) + 0.01;
    requiredY(1,1) = requiredY(1,1) - 0.01;
    requiredY(2,1) = requiredY(2,1) + 0.01;
    results_file = [results_file, 'Pr'];
else
    opticalCavity = [requiredX(:,1) requiredY(:,1)];
end
%% Create coordinates
coordinatesX = transpose(linspace(xRange(1), xRange(2), NX));
coordinatesX = sort(unique([coordinatesX; addedX(:); requiredX(isfinite(requiredX)); opticalCavity(:,1)]));
coordinatesY = transpose(linspace(yRange(1), yRange(2), NY));
coordinatesY = sort(unique([coordinatesY; addedY(:); requiredY(isfinite(requiredY)); opticalCavity(:,2)]));
NX = length(coordinatesX);
NY = length(coordinatesY);
coordinates = [repmat(coordinatesX, NY, 1), reshape(transpose(repmat(coordinatesY, 1, NX)), [], 1)];
x = coordinates(:,1);
y = coordinates(:,2);
%% Create elements
elementsBlock = [1 2 NX+2; 2 3 NX+2; 1 NX+2 NX+1; 3 NX+3 NX+2; NX+1 NX+2 2*NX+1; NX+2 NX+3 2*NX+3; NX+2 2*NX+2 2*NX+1; NX+2 2*NX+3 2*NX+2];
e2p = nan(2*(NX-1)*(NY-1),3);
if mod(NX, 2) ~= 1 || mod(NY, 2) ~= 1
    error('NX and NY must be odd');
end
ijCounter = 0;
for i=1:(NX-1)/2
    for j=1:(NY-1)/2
        ijCounter = ijCounter + 1;
        e2p(8*(ijCounter-1)+1:8*ijCounter,:) = elementsBlock+(i-1)*2+(j-1)*NX*2;
    end
end



[TriInfo, Transformation] = MeshCreateMatrices1(x,y,e2p,requiredX,requiredY,sizePhi,opticalCavity);
[Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
matrices.H1scalProlong =  TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
save('Data', 'Transformation', 'TriInfo', 'matrices');
mkdir(results_file);
movefile('Data.mat',results_file);
%% plot the data
if plotResults
    figures = MeshCreateFigures(TriInfo);
    for i=1:length(figures)
        figureName = ['Visualization', int2str(i), '.jpg'];
        saveas(figures{i}, figureName);
        movefile(figureName, results_file);
    end
end
