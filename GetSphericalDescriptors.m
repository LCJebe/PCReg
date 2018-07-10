clear all
close all

%% TODOs
% - reject spherical features (DONE: PCA variances) 
% - use alignment strategy  (DONE: PCA w/ sign disambiguition)
% - switch to spherical harmonics

ALIGN_POINTS = true;

%% define method for sampling points (keypoint detection)
% 'UNIFORM' for regular grid (uniform sampling)
% 'RANDOM' for random points in point cloud
% 'ALL' for all points in point cloud
method = 'UNIFORM';

%% READ in aligned surface and model crop
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'Surface_DS_aligned.pcd'));
pcModel = pcread(strcat(path, 'GoodCrop.pcd'));
pcRand = pcread(strcat(path, 'RandCrop.pcd'));

%% define locations of spheres for local descriptors
if strcmp(method, 'UNIFORM')
    XLimits = [min(pcSurface.XLimits(1), pcModel.XLimits(1)), ...
                max(pcSurface.XLimits(2), pcModel.XLimits(2))];
    YLimits = [min(pcSurface.YLimits(1), pcModel.YLimits(1)), ...
                max(pcSurface.YLimits(2), pcModel.YLimits(2))];
    ZLimits = [min(pcSurface.ZLimits(1), pcModel.ZLimits(1)), ...
                max(pcSurface.ZLimits(2), pcModel.ZLimits(2))];

    % create meshgrid of points
    d = 0.4; % spacing of spheres along each axis
    x = XLimits(1):d:XLimits(2);
    y = YLimits(1):d:YLimits(2);
    z = ZLimits(1):d:ZLimits(2);
    [X,Y,Z] = meshgrid(x,y,z);
    sample_pts = cat(4, X, Y, Z);
    sample_pts = reshape(sample_pts, [], 3);
end

%% get local descriptors for each model

% specify minimum number of points that has to be in sphere
min_pts = 50;
% specify radius of sphere
R = 1.5;
% specify the reject threshold for eccentricity (covar-eigenvalues), value
% must be >= 1 
thVar = [1.5, 1.5];

% specify number of nearest neighbors (KNN) to use for local reference
% frame. Number should be <= min_points, or write 'all'
% if k is 'all', then points need not be sorted - faster. 
k = 'all';

pts = pcModel.Location;
[featModel, descModel] = getMomentDescriptors(pts, sample_pts, min_pts, R, thVar, ALIGN_POINTS, k);

pts = pcRand.Location;
[featRand, descRand] = getMomentDescriptors(pts, sample_pts, min_pts, R, thVar, ALIGN_POINTS, k);

pts = pcSurface.Location;
[featSurface, descSurface] = getMomentDescriptors(pts, sample_pts, min_pts, R, thVar, ALIGN_POINTS, k);

%% Apply weights to descriptors
% first step: weight descriptors
weights.M1 = 0.2; % 0.2
weights.M2 = 1; % 1
weights.M3 = 1.8; % 1.8

descSurfaceW = applyWeight(descSurface, weights);
descModelW = applyWeight(descModel, weights);
descRandW = applyWeight(descRand, weights);
%% Match features between Surface and Model / Random Crop

% Define matching algorithm parameters
par.Method = 'Exhaustive'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold =  0.3; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.95; % 0.6 (default) nearest neighbor ambiguity rejection
par.Metric =  'SSD'; % SSD (default) for L2, SAD for L1
par.Unique = true; % true: 1-to-1 mapping only, else set false (default)



matchesModel = matchFeatures(descSurfaceW, descModelW, ...
        'Method', par.Method, ...
        'MatchThreshold', par.MatchThreshold, ... 
        'MaxRatio', par.MaxRatio, ... 
        'Metric', par.Metric, ...
        'Unique', par.Unique); 

matchesRand = matchFeatures(descSurfaceW, descRandW, ...
        'Method', par.Method, ...
        'MatchThreshold', par.MatchThreshold, ... 
        'MaxRatio', par.MaxRatio, ... 
        'Metric', par.Metric, ...
        'Unique', par.Unique); 
    
%% Get distance between matching points
% this makes sense, because the pointclouds are already aligned. The
% distance between matching INLIERS of the model will thus be small and a
% simple distance threshold can be used to determine whether a match is an
% inlier 

% maxDist specifies matching distance to count inliers 
maxDist = 1; 

% matches of surface and model
loc1M = featModel(matchesModel(:, 2), :);
loc1S = featSurface(matchesModel(:, 1), :);
d1 = vecnorm(loc1M - loc1S, 2, 2);
inliers1 = length(find(d1 < maxDist));

% precision
inliersPrecision = inliers1/size(d1, 1)*100;

% matches of surface and random crop
loc2M = featRand(matchesRand(:, 2), :);
loc2S = featSurface(matchesRand(:, 1), :);
d2 = vecnorm(loc2M - loc2S, 2, 2);
inliers2 = length(find(d2 <= maxDist));

                   

%% Function: get points within sphere
% returns only the points from pts that are within a certain radius R
% of center point c if there are at least min_points in it
% RETURNS: Points RELATIVE to c, optinally sorted by distance. 
function pts_sphere = getLocalPoints(pts, R, c, min_points, SORT_POINTS)

    % first: quickly crop cube around the center
    xLim = c(1) + [-R, R];
    yLim = c(2) + [-R, R];
    zLim = c(3) + [-R, R];
    mask = pts(:, 1) > xLim(1) & pts(:, 1) < xLim(2) ...
         & pts(:, 2) > yLim(1) & pts(:, 2) < yLim(2) ...
         & pts(:, 3) > zLim(1) & pts(:, 3) < zLim(2);
    mask = cat(2, mask, mask, mask);
    pts_cube = reshape(pts(mask), [], 3);    
    
    if size(pts_cube, 1) < min_points
        pts_sphere = [];
    else
        % then: get distances for each remaining point and only keep the points
        % that are within the radius
        pts_rel = pts_cube - c;
        dists = vecnorm(pts_rel, 2, 2);
        mask = dists < R;
        dists_remaining = dists(find(mask));
        mask = cat(2, mask, mask, mask);
        pts_sphere = reshape(pts_rel(mask), [], 3);    
        
        % reject if not enough points found
        if size(pts_sphere, 1) < min_points
            pts_sphere = [];
        elseif SORT_POINTS % order points by distance
            [~, I] = sort(dists_remaining);
            pts_sphere = pts_sphere(I, :);
        end
    end
    
end

%% helpter function that applies weight to descriptor
function descW = applyWeight(desc, weights)
    descW = desc;
    descW(:, 1:6) = desc(:, 1:6)*weights.M1;
    descW(:, 7:12) = desc(:, 7:12)*weights.M2;
    descW(:, 13:22) = desc(:, 13:22)*weights.M3;
end