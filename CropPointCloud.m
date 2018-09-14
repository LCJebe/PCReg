clear all
close all

%% DESCRIPTION
% this script crops out part of the pointcloud as a sphere. 
% this makes sense, assuming we don't know the true alignment. 
% We need to specify the center and the crop radius. The radius is flexible,
% and a parameter for the experiement. To include all possibilities, select 
% R = len(diag(cuboid)) + 3.5 (descriptor radius) ("ideal")

%% options
center = [52, 12, 44]; % specify center for crop
cropsize = 'smaller'; % 'smaller', 'smallest', or 'ideal'
name = 'NewBadCrop3.pcd'; % name, including ".pcd"


%% read pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothColorUp3.pcd'));

% load limits of aligned surface for crop-reference!
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));
xLimS = pcSurface.XLimits;
yLimS = pcSurface.YLimits;
zLimS = pcSurface.ZLimits;

% get length of cuboid diagonal
diagSurface = norm([xLimS(2) - xLimS(1), ...
                    yLimS(2) - yLimS(1), ...
                    zLimS(2) - zLimS(1)]);

%% specify final crop radius
if strcmp(cropsize, 'ideal')
    R_crop = diagSurface/2 + 3.5;
elseif strcmp(cropsize, 'smaller')
    R_crop = diagSurface/2;
elseif strcmp(cropsize, 'smallest')
    R_crop = diagSurface/2 - 3.5;
else
    error('please choose a valid cropsize option')
end

%% create spherical crop of pointcloud
pts = pcModel.Location;
colors = pcModel.Color;

pts_rel = pts - center;
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
pcwrite(pcCrop, strcat('Data/PointClouds/', name));