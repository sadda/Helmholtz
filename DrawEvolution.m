workingDir = fullfile('Results', 'Results_Alternating2_RegThetaL1=0.600', 'Pictures');
imageNames = dir(fullfile(workingDir, '*.jpg'));
imageNames = {imageNames.name}';
imageNames = imageNames(1:length(imageNames)/2);





imageNames2 = imageNames;
for ii=10:17
    imageNames{ii} = imageNames2{ii+9};
end
for ii=18:26
    imageNames{ii} = imageNames2{ii-8};
end    
    
    
% imageNames




outputVideo = VideoWriter('Evolution.avi');
outputVideo.FrameRate = 3;
open(outputVideo);

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir, imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)