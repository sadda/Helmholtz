clear all;
close all;

% Always check the mesh visually before working with it :)

fixGermanium = 0;
plotResults  = 1;
xRange       = [-2 2];
yRange       = [0 3];
% NX           = 33;
% NY           = 25;
NX = 17;
NY = 13;
sizePhi      = 4;
if fixGermanium
    requiredX    = [-0.75 -0.75 -2 -2; 0.75 0.75 2 2];
    requiredY    = [1.1 1.48 0 2.5; 1.4 1.8 0.96 3];
    results_file = 'MeshesCreated/Mesh_GeFixed';
    addedX       = [];
    addedY       = [1.445];
else
    requiredX    = [Inf Inf -2 -2; -Inf -Inf 2 2];
    % requiredY    = [Inf Inf 0 2.5; -Inf -Inf 3/24*7 3];
    requiredY    = [Inf Inf 0 2.5; -Inf -Inf 1 3];
    results_file = 'MeshesCreated/Mesh_GeFree';
    addedX       = [];
    addedY       = [];
end

%% Create coordinates
coordinatesX = transpose(linspace(xRange(1), xRange(2), NX));
coordinatesX = sort(unique([coordinatesX; addedX(:); requiredX(isfinite(requiredX))]));
coordinatesY = transpose(linspace(yRange(1), yRange(2), NY));
coordinatesY = sort(unique([coordinatesY; addedY(:); requiredY(isfinite(requiredY))]));
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
% Only necessary for the refining and coarsening functions which I use
i1        = 1+8*(0:(size(e2p,1)/8-1));
i2        = 2+8*(0:(size(e2p,1)/8-1));
i4        = 4+8*(0:(size(e2p,1)/8-1));
i5        = 5+8*(0:(size(e2p,1)/8-1));
i6        = 6+8*(0:(size(e2p,1)/8-1));
i7        = 7+8*(0:(size(e2p,1)/8-1));
e2p(i1,:) = e2p(i1,[3 1 2]);
e2p(i2,:) = e2p(i2,[2 3 1]);
e2p(i4,:) = e2p(i4,[3 1 2]);
e2p(i5,:) = e2p(i5,[2 3 1]);
e2p(i6,:) = e2p(i6,[3 1 2]);
e2p(i7,:) = e2p(i7,[3 1 2]);

[TriInfo, Transformation] = MeshCreateMatrices1(x,y,e2p,requiredX,requiredY,sizePhi);
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
