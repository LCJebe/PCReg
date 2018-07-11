clear all
close all

%% TODOs
% - reject spherical features (DONE: PCA variances) 
% - use alignment strategy  (DONE: PCA w/ sign disambiguition)

%% define method for sampling points (keypoint detection)
% 'UNIFORM' for regular grid (uniform sampling)
% 'RANDOM' for random points in point cloud
% 'ALL' for all points in point cloud
sampling_method = 'UNIFORM';

%% READ in aligned surface and model crop
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'Surface_DS_aligned.pcd'));
pcModel = pcread(strcat(path, 'GoodCrop.pcd'));
pcRand = pcread(strcat(path, 'RandCrop.pcd'));

%% define locations of spheres for local descriptors
if strcmp(sampling_method, 'UNIFORM')
    XLimits = [min(pcSurface.XLimits(1), pcModel.XLimits(1)), ...
                max(pcSurface.XLimits(2), pcModel.XLimits(2))];
    YLimits = [min(pcSurface.YLimits(1), pcModel.YLimits(1)), ...
                max(pcSurface.YLimits(2), pcModel.YLimits(2))];
    ZLimits = [min(pcSurface.ZLimits(1), pcModel.ZLimits(1)), ...
                max(pcSurface.ZLimits(2), pcModel.ZLimits(2))];

    % create meshgrid of points
    d = 0.3; % 0.4, spacing of spheres along each axis
    x = XLimits(1):d:XLimits(2);
    y = YLimits(1):d:YLimits(2);
    z = ZLimits(1):d:ZLimits(2);
    [X,Y,Z] = meshgrid(x,y,z);
    sample_pts = cat(4, X, Y, Z);
    sample_pts = reshape(sample_pts, [], 3);
end

%% get local descriptors for each model

% define whether to use 'Moment'- or 'Rotational' Descriptor
descriptor = 'Moment';

% if descriptor is 'Moment' optionally align to a local reference frame
ALIGN_POINTS = true;

% specify minimum number of points that has to be in sphere
min_pts = 50; % 50
% specify radius of sphere
R = 1.5; % 1.5
% specify the reject threshold for eccentricity (covar-eigenvalues), value
% must be >= 1 
thVar = [1.5, 1.5];

% specify number of nearest neighbors (KNN) to use for local reference
% frame. Number should be <= min_points, or write 'all'
% if k is 'all', then points need not be sorted - faster. 
k = 'all';

% get points from pointclouds
ptsModel = pcModel.Location;
ptsRand = pcRand.Location;
ptsSurface = pcSurface.Location;

% calculate the descriptos with the specified method
if strcmp(descriptor, 'Moment')    
    [featModel, descModel] = getMomentDescriptors(ptsModel, sample_pts, min_pts, R, thVar, ALIGN_POINTS, k);
    [featRand, descRand] = getMomentDescriptors(ptsRand, sample_pts, min_pts, R, thVar, ALIGN_POINTS, k);
    [featSurface, descSurface] = getMomentDescriptors(ptsSurface, sample_pts, min_pts, R, thVar, ALIGN_POINTS, k);        
elseif strcmp(descriptor, 'Rotational')    
    [featModel, descModel] = getRotationalDescriptors(ptsModel, sample_pts, min_pts, R);
    [featRand, descRand] = getRotationalDescriptors(ptsRand, sample_pts, min_pts, R);
    [featSurface, descSurface] = getRotationalDescriptors(ptsSurface, sample_pts, min_pts, R);
else
    error('Descriptor method must be either Moment or Rotational');
end

    
%% Apply weights to descriptors
if strcmp(descriptor, 'Moment')    
    weights.M1 = 0; % 0.2
    weights.M2 = 0.3; % 1
    weights.M3 = 0.8; % 1.8
    weights.M4 = 0.5; % 0.5

    descSurfaceW = applyMomentWeight(descSurface, weights);
    descModelW = applyMomentWeight(descModel, weights);
    descRandW = applyMomentWeight(descRand, weights);
else
    descSurfaceW = descSurface;
    descModelW = descModel;
    descRandW = descRand;
end
    

%% Match features between Surface and Model / Random Crop

% Define matching algorithm parameters
par.Method = 'Exhaustive'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold =  1; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.8; % 0.6 (default) nearest neighbor ambiguity rejection
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

%% helper function that applies weight to descriptor
function descW = applyMomentWeight(desc, weights)
    descW = desc;
    descW(:, 1:6) = desc(:, 1:6)*weights.M1;
    descW(:, 7:12) = desc(:, 7:12)*weights.M2;
    descW(:, 13:22) = desc(:, 13:22)*weights.M3;
    descW(:, 23:31) = desc(:, 23:31)*weights.M4;
end