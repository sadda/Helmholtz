function CropImage(imageName, imageName2)
    
    if nargin < 2 || ~ischar(imageName2)
        [filePart1, filePart2, filePart3] = fileparts(imageName);
        imageName2 = [fullfile(filePart1, filePart2), '_Cropped', filePart3];
    end
    
    image = imread(imageName);
    
    imageMod = image;
    if max(imageMod(:)) <= 1
        imageMod = 255*imageMod;
    end
    if ismatrix(imageMod)
        imageMod = repmat(imageMod, 1, 1, 3);
    end
    imageSum = sum(imageMod, 3);
    nonzero1 = find(sum(imageSum, 2) ~= size(imageSum,2)*255*3);
    nonzero2 = find(sum(imageSum, 1) ~= size(imageSum,1)*255*3);
    if ismatrix(image)
        image    = image(nonzero1(1):nonzero1(end), nonzero2(1):nonzero2(end));
    else
        image    = image(nonzero1(1):nonzero1(end), nonzero2(1):nonzero2(end), :);
    end
    
    imwrite(image, imageName2);
end