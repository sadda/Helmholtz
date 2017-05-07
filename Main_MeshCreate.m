% Creates mesh NX*XY nodes and the corresponding matrices.

clear all;
close all;

xRange       = [-2 2];
yRange       = [0 3];
NX           = 33;
NY           = 25;
sizePhi      = 4;
plotResults  = 1;

requiredX    = [-0.125 -0.75 -2 -2; 0.125 0.75 2 2]; % x range for box Pi_i
requiredY    = [1 1.5 0 2.5; 1.49 1.75 0.99 3]; % y range for box Pi_i
results_file = 'MeshesCreated/Mesh_Final';
addedX       = [];
addedY       = [];

%% Create coordinates

coordinatesX = transpose(linspace(xRange(1), xRange(2), NX));
coordinatesX = sort(unique([coordinatesX; addedX(:)]));
coordinatesY = transpose(linspace(yRange(1), yRange(2), NY));
coordinatesY = sort(unique([coordinatesY; addedY(:)]));
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
% For the refinement and coarsening from P1AFEM, some of the elements need to be rotated
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

%% Create matrices and save the data

[TriInfo, Transformation]           = MeshCreateMatrices1(x,y,e2p,requiredX,requiredY,sizePhi);
[Transformation, TriInfo, matrices] = MeshCreateMatrices2(Transformation, TriInfo);
matrices.H1scalProlong              =  TriInfo.phiProlongationMatrix6'*matrices.H1scal*TriInfo.phiProlongationMatrix6;
save('Data', 'Transformation', 'TriInfo', 'matrices');
mkdir(results_file);
movefile('Data.mat',results_file);

%% Plot the data

if plotResults
    figures = MeshCreateFigures(TriInfo);
    for i=1:length(figures)
        figureName = ['Visualization', int2str(i), '.jpg'];
        saveas(figures{i}, figureName);
        movefile(figureName, results_file);
    end
end
