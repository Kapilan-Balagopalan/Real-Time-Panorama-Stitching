fileNames={'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\a.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\b.jpg',...
    'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\c.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\d.jpg',...
    'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\e.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\6.jpg',...
    'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\7.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\8.jpg',...
    'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\9.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\10.jpg',...
    'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\11.jpg'
    };
% vid = videoinput('winvideo', 1, 'MJPG_640x480');
% src = getselectedsource(vid);
% 
% 
%  
%  vid.FramesPerTrigger = 1600;
% preview(vid);
% sad=input('tell me');
% start(vid);
% % pause(10);
% 
% 
% data = getdata(vid,120);
% size(data)
% stop(vid);
% delete(vid);
% %image = imread('Clus4.jpg');
% for i=1:5
%     
% image=(data(:,:,:,24*i));
% image = imresize(image,1);
% imwrite(image,fileNames{i});
% end
% % Load images.
% 
% url = 'http://192.168.43.1:8080/photo.jpg';
% ss  = imread(url);
% 
% %fh = image(ss);
% q=1;
% while(q<6)
%     ss  = imread(url);
% %     imshow(ss);
%      imwrite(ss,fileNames{q});
%     pause(1);
%     q=q+1;
%     
% end
% 
% disp('taken');
% buildingDir = fullfile(toolboxdir('vision'), 'visiondata', 'building');
% dirOutput=dir(fullfile(buildingDir,'*.png'));
% %buildingScene = imageDatastore(buildingDir);
% fileNames = {dirOutput.name}';
% imshow('C:\Program Files\MATLAB\R2013a\toolbox\vision\visiondata\building\1.png')
% fileNames={'C:\Program Files\MATLAB\R2013a\toolbox\vision\visiondata\building\1.png','C:\Program Files\MATLAB\R2013a\toolbox\vision\visiondata\building\2.png',...
%     'C:\Program Files\MATLAB\R2013a\toolbox\vision\visiondata\building\3.png','C:\Program Files\MATLAB\R2013a\toolbox\vision\visiondata\building\4.png',...
%     'C:\Program Files\MATLAB\R2013a\toolbox\vision\visiondata\building\5.png'}
% % 
% fileNames={'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\1.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\2.jpg',...
%     'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\3.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\4.jpg',...
%     'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\5.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\6.jpg',...
%     'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\7.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\8.jpg',...
%     'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\9.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\10.jpg',...
%     'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\11.jpg'
%     };
% 
% fileNames={'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\trees_007.jpg',...
%     'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\trees_008.jpg'
%     };


%  fileNames={'C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\CMU_left.jpg','C:\Users\HP\Desktop\Panorama-Mosaicing-master\Photos\CMU_right.jpg'...
%    };
% Display images to be stitched
%montage(fileNames);

% Read the first image from the image set.
%I = readimage(buildingScene, 1);
I=imread(fileNames{1});
 % I=imresize(I,0.125);
%imshow(I);
% Initialize features for I(1)
grayImage = rgb2gray(I);
points = detectSURFFeatures(grayImage);
[features, points] = extractFeatures(grayImage, points);
% imshow(I); hold on;
% plot(points.selectStrongest(20));

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = length(fileNames);

numImages =4;
tforms(numImages) = projective2d(eye(3));

% Iterate over remaining image pairs
for n = 2:numImages 

    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;

    % Read I(n).
   % I = readimage(buildingScene, n)
    I=imread(fileNames{n});
    I1=I;
    Ip=imread(fileNames{n-1});
     % I=imresize(I,0.125);
%       Ip=imresize(Ip,0.5);
    % Detect and extract SURF features for I(n).
    grayImage = rgb2gray(I);
    points = detectSURFFeatures(grayImage);
    [features, points] = extractFeatures(grayImage, points);

    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious,'Prenormalized',false);

    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);
%     showMatchedFeatures(I1,Ip,matchedPoints,matchedPointsPrev);
   
    % Estimate the transformation between I(n) and I(n-1).
    tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);

    % Compute T(1) * ... * T(n-1) * T(n)
    tforms(n).T = tforms(n-1).T * tforms(n).T;
end


imageSize = size(I);  % all the images are the same size

% Compute the output limits  for each transform
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end

avgXLim = mean(xlim, 2);

[~, idx] = sort(avgXLim);

centerIdx = floor((numel(tforms)+1)/2);

centerImageIdx = idx(centerIdx);

Tinv = invert(tforms(centerImageIdx));

for i = 1:numel(tforms)
    tforms(i).T = Tinv.T * tforms(i).T;
end

for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
end
Ycoord=zeros(numel(tforms),4,2);

for i = 1:numel(tforms)
     cornerCoord(i,:,:)=transformPointsForward(tforms(i),[0,0;0,imageSize(2);imageSize(1),0;imageSize(1),imageSize(2)]);

end

[valymax idxymax] =max(max(cornerCoord(:,:,2).'));

[valymin idxymin] =min(min(cornerCoord(:,:,2).'));

[valxmax idxxmax] =max(max(cornerCoord(:,:,1).'));

[valxmin idxxmin] =min(min(cornerCoord(:,:,1).'));


deltay1=abs(cornerCoord(idxymin,1,2)-cornerCoord(idxymin,3,2));
deltay2=abs(cornerCoord(idxymax,2,2)-cornerCoord(idxymax,4,2));

deltax1=abs(cornerCoord(idxxmin,1,1)-cornerCoord(idxxmin,2,1));
deltax2=abs(cornerCoord(idxxmax,3,1)-cornerCoord(idxxmax,4,1));
% Find the minimum and maximum output limits
xMin = min([1; xlim(:)]);
xMax = max([imageSize(2); xlim(:)]);

yMin = min([1; ylim(:)]);
yMax = max([imageSize(1); ylim(:)]);

cutx1=round(max(xlim(:,1)));
cutx2=round(min(xlim(:,2)));

cuty1=round(max(ylim(:,1)));
cuty2=round(max(ylim(:,2)));
% Width and height of panorama.
width  = round(xMax - xMin);
height = round(yMax - yMin);

% Initialize the "empty" panorama.
panorama = zeros([height width 3], 'like', I);


blender = vision.AlphaBlender('Operation','Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
xLimits = [xMin xMax];
yLimits = [yMin yMax];
panoramaView = imref2d([height width], xLimits, yLimits);

%  warpedImage = imwarp(I, tforms(1), 'OutputView', panoramaView);
%  size(warpedImage)
%  imshow(warpedImage);
%  return
% Create the panorama.
setSize=imageSize;
e1=15;
for i = 1:numImages

    %I = readimage(buildingScene, i);
     I=imread(fileNames{i});
     % I=imresize(I,0.125);
    % Transform I into the panorama.
    warpedImage = imwarp(I, tforms(i), 'OutputView', panoramaView);
    
    setSize=[setSize; size(warpedImage)];
    % Generate a binary mask.
    mask = imwarp(true(size(I,1),size(I,2)), tforms(i), 'OutputView', panoramaView);
    
    
    % Overlay the warpedImage onto the panorama.
    panorama = step(blender, panorama, warpedImage,mask);
end
imwrite(panorama,'panaorama_raw.jpg');
% figure,
% imshow(panorama)
panorama=panorama(deltay1+e1:end-deltay2-e1,deltax1+e1:end-deltax2-e1,:);
figure
imshow(panorama)
imwrite(panorama,'panaorama.jpg');
