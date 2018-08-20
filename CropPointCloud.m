clear all
close all

%% read pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothColor.pcd'));

% load limits of aligned surface for crop-reference!
pcSurface = pcread(strcat(path, 'Surface_DS4_alignedM.pcd'));
xLimS = pcSurface.XLimits;
yLimS = pcSurface.YLimits;
zLimS = pcSurface.ZLimits;

% optional: show pointcloud
% figure()
% pcshow(pcModel, 'MarkerSize', 50);
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

%% specify limits in [mm]
margin = 3.5;
xLim = xLimS + [-margin, margin];
yLim = yLimS + [-margin, margin];
zLim = zLimS + [-margin, margin];

%% create crop of pointcloud
pts = pcModel.Location;
colors = pcModel.Color;
mask = pts(:, 1) > xLim(1) & pts(:, 1) < xLim(2) ...
     & pts(:, 2) > yLim(1) & pts(:, 2) < yLim(2) ...
     & pts(:, 3) > zLim(1) & pts(:, 3) < zLim(2);
mask = cat(2, mask, mask, mask);
pts_crop = reshape(pts(mask), [], 3);  
colors_crop = reshape(colors(mask), [], 3);  

pcCrop = pointCloud(pts_crop, 'Color', colors_crop);

% optional: show pointcloud
% figure()
% pcshow(pcCrop, 'MarkerSize', 50);
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

%% save the crop
pcwrite(pcCrop, 'Data/PointClouds/GoodCropSmooth_large.pcd');