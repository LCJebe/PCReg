clear all
close all

%% read pointcloud
path = 'Data/PointClouds/Model.pcd';
pcModel = pcread(path);

% optional: show pointcloud
figure()
pcshow(pcModel, 'MarkerSize', 50);
xlabel('X');
ylabel('Y');
zlabel('Z');

%% specify limits in [mm]
xLim = [30, 50];
yLim = [15, 35];
zLim = [15, 35];

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
figure()
pcshow(pcCrop, 'MarkerSize', 50);
xlabel('X');
ylabel('Y');
zlabel('Z');