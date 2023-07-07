% cropPlot automatically crops the whitespace out of any image provided
% Needs the full path to the image. Automatically saves the cropped image
% with an added extension (-crop)
function cropPlot(imPath)
    im = imread(imPath);
    
    [~, ~, ext] = fileparts(imPath);
    % Get mean value across 3-dimension (255 is white)
    meanIm = mean(im,3);
    
    % Look at the columns
    meanCol  = mean(meanIm,2);
    getWhite = meanCol ~= 255;
    getCol   = [find(getWhite, 1) , find(getWhite, 1, 'last')];
    
    % look at the rows
    meanRow  = mean(meanIm,1);
    getWhite = meanRow ~= 255;
    getRow   = [find(getWhite, 1) , find(getWhite, 1, 'last')];
    
    % Get the size of the image
    imSize = [getCol(2) - getCol(1), getRow(2) - getRow(1)];
    
    % Crop the image
    cropIm = imcrop(im,[getRow(1), getCol(1), imSize(2), imSize(1)]);
    
    % Save the plot
    imwrite(cropIm, replace(imPath, ext, ['-crop' ext ]))