clear all
close all


%% DESCRIPTION
% this script crops out the good crop (optionally bad crop) as a sphere. 
% this makes sense, assuming we don't know the true alignment. 
% We need to specify the center and the crop radius. For the ideal good
% crop, the center has to be based on the center of the cuboid (Quader) of
% the aligned Surface. The radius is flexible, and a parameter for the
% experiement. To include all possibilities, select 
% R = len(diag(cuboid)) + 3.5 (descriptor radius)

% options
RAND_CROP = true;

%% read pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothColorUp3.pcd'));

% load limits of aligned surface for crop-reference!
pcSurface = pcread(strcat(path, 'Surface_DS3_alignedM.pcd'));
xLimS = pcSurface.XLimits;
yLimS = pcSurface.YLimits;
zLimS = pcSurface.ZLimits;

% get center of aligned Surface
centerSurface = [(xLimS(2) + xLimS(1))/2, ...
                 (yLimS(2) + yLimS(1))/2, ...
                 (zLimS(2) + zLimS(1))/2];
             
% get length of cuboid diagonal
diagSurface = norm([xLimS(2) - xLimS(1), ...
                    yLimS(2) - yLimS(1), ...
                    zLimS(2) - zLimS(1)]);

% optional: show pointcloud
% figure()
% pcshow(pcModel, 'MarkerSize', 50);
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

%% specify final crop radius
%R_crop = diagSurface/2 + 3.5; % diagSurface/2 + 3.5 for "ideal" crop
%(16.5mm)
%R_crop = diagSurface/2; % smaller
R_crop = diagSurface/2 - 3.5; % smallest


% for bad crop only: shift point cloud by xShift TODOOOOOO!! (unfinished)
margin = 3.5;
if RAND_CROP
    xShift = diagSurface/2 + R_crop;
    T = eye(4);
    T(4, 1) = xShift;
    pcModel = pctransform(pcModel, affine3d(T));
end

%% create spherical crop of pointcloud
pts = pcModel.Location;
colors = pcModel.Color;

pts_rel = pts - centerSurface;
dists = vecnorm(pts_rel, 2, 2);
mask = dists < R_crop;
mask = cat(2, mask, mask, mask);
pts_crop = reshape(pts(mask), [], 3);  
colors_crop = reshape(colors(mask), [], 3);  

pcCrop = pointCloud(pts_crop, 'Color', colors_crop);

% optional: show pointcloud
figure()
pcshow(pcCrop, 'MarkerSize', 50);
xlabel('X');
ylabel('Y');
zlabel('Z');

%% save the crop
pcwrite(pcCrop, 'Data/PointClouds/RandCropSpherical_smallest.pcd');