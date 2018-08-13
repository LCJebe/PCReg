clear all
close all

%% TODOs
% - implement keypoint detector

%% define method for sampling points (keypoint detection)
% 'UNIFORM' for regular grid (uniform sampling)
% 'UNIFORM_SAME' for the same regular grid on all point clouds
% 'RANDOM_POINTS' for random points in point cloud
% 'RANDOM_UNIFORM' for random uniformly distributed points
% 'ALL' for all points in point cloud
% 'SELECTIVE' for certain distinct points in point cloud (not implemented
% yet)
sampling_method = 'RANDOM_UNIFORM';

% if 'UNIFORM' or 'UNIFORM_SAME' or 'RANDOM_UNIFORM' specify density of sphere-grid
d = 0.2; % 0.4, spacing of spheres along each axis (0.3 is better, 0.2 takes forever)

%if 'RANDOM_POINTS' specify fraction of points considered
sample_frac = 0.2;

%% READ in aligned surface and model crop
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'Surface_DS2_aligned.pcd'));
pcModel = pcread(strcat(path, 'GoodCrop_v2.pcd'));
pcRand = pcread(strcat(path, 'RandCrop_v2.pcd'));

% recommended: center point clouds
SHIFT = false; % 'true' destroys alignment, if aligned

% turn off to omit matching with random crop
RAND_CROP = true;

if SHIFT
    pcSurface = centerPointCloud(pcSurface);
    pcModel = centerPointCloud(pcModel);
    pcRand = centerPointCloud(pcRand);
end

% DEBUG: manually align now REMOVE for normal experiments!!!
%load('Tform_align2.mat');
%pcSurface = pctransform(pcSurface, affine3d(Tform_align2));
%pcSurface = centerPointCloud(pcSurface);

%% define locations of spheres for local descriptors
if strcmp(sampling_method, 'UNIFORM_SAME')

    XLimits = [min([pcSurface.XLimits(1), pcModel.XLimits(1), pcRand.XLimits(1)]), ...
                max([pcSurface.XLimits(2), pcModel.XLimits(2), pcRand.XLimits(2)])];
    YLimits = [min([pcSurface.YLimits(1), pcModel.YLimits(1), pcRand.YLimits(1)]), ...
                max([pcSurface.YLimits(2), pcModel.YLimits(2), pcRand.YLimits(2)])];
    ZLimits = [min([pcSurface.ZLimits(1), pcModel.ZLimits(1), pcRand.ZLimits(1)]), ...
                max([pcSurface.ZLimits(2), pcModel.ZLimits(2), pcRand.ZLimits(2)])];

    % create meshgrid of points    
    x = XLimits(1):d:XLimits(2);
    y = YLimits(1):d:YLimits(2);
    z = ZLimits(1):d:ZLimits(2);
    [X,Y,Z] = meshgrid(x,y,z);
    sample_pts = cat(4, X, Y, Z);
    sample_pts = reshape(sample_pts, [], 3);
    
    sample_ptsSurface = sample_pts;
    sample_ptsModel = sample_pts;
    sample_ptsRand = sample_pts;
    
elseif strcmp(sampling_method, 'UNIFORM')
    
    sample_ptsSurface = pcUniformSamples(pcSurface, d);
    sample_ptsModel = pcUniformSamples(pcModel, d);
    sample_ptsRand = pcUniformSamples(pcRand, d);
    
elseif strcmp(sampling_method, 'RANDOM_UNIFORM')
    
   sample_ptsSurface = pcRandomUniformSamples(pcSurface, d);
   sample_ptsModel = pcRandomUniformSamples(pcModel, d);
   sample_ptsRand = pcRandomUniformSamples(pcRand, d);
    
elseif strcmp(sampling_method, 'RANDOM_POINTS')
    % Surface
    nPoints = size(pcSurface.Location, 1);
    sIdx = randsample(nPoints,floor(nPoints*sample_frac));
    sample_ptsSurface = pcSurface.Location(sIdx, :);
    
    % Model
    nPoints = size(pcModel.Location, 1);
    sIdx = randsample(nPoints,floor(nPoints*sample_frac));
    sample_ptsModel = pcModel.Location(sIdx, :);
    
    % Rand
    nPoints = size(pcRand.Location, 1);
    sIdx = randsample(nPoints,floor(nPoints*sample_frac));
    sample_ptsRand = pcRand.Location(sIdx, :);
    
elseif strcmp(sampling_method, 'ALL')
    sample_ptsSurface = pcSurface.Location;
    sample_ptsModel = pcModel.Location;
    sample_ptsRand = pcRand.Location;
else
    error('Specified Method "%s" not yet implemented or does not exist', sampling_method);
end

%% get local descriptors for each model

% define which descriptor to use
% 'Moment' --> using 1st to 4th order moments with LRF
% 'Rotational' --> inherently rotation invariant (radial frequencies)
% 'Histogram' --> Spacial spherical histogram
descriptor = 'Moment';

% if descriptor is 'Moment' optionally align to a local reference frame
ALIGN_POINTS = true;

% optional: don't use first moment and center by this instead
CENTER = false;

% specify minimum number of points that has to be in sphere
min_pts = 101; % 51 / 101

% specify radius of spheres (local descriptor neighborhood)
R = 2.5; % 1.5 / 2.5

% specify the reject threshold for eccentricity (covar-eigenvalues), value
% must be >= 1 
thVar = [3, 1]; % [1.5, 1.5] 

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
    [featModel, descModel, angModel] = getMomentDescriptors(ptsModel, sample_ptsModel, min_pts, R, thVar, ALIGN_POINTS, CENTER, k);
    if RAND_CROP
        [featRand, descRand, angRand] = getMomentDescriptors(ptsRand, sample_ptsRand, min_pts, R, thVar, ALIGN_POINTS, CENTER, k);
    end
    [featSurface, descSurface, angSurface] = getMomentDescriptors(ptsSurface, sample_ptsSurface, min_pts, R, thVar, ALIGN_POINTS, CENTER, k);        
elseif strcmp(descriptor, 'Rotational')    
    [featModel, descModel] = getRotationalDescriptors(ptsModel, sample_ptsModel, min_pts, R);
    if RAND_CROP
        [featRand, descRand] = getRotationalDescriptors(ptsRand, sample_ptsRand, min_pts, R);
    end
    [featSurface, descSurface] = getRotationalDescriptors(ptsSurface, sample_ptsSurface, min_pts, R);
elseif strcmp(descriptor, 'Histogram')    
    [featModel, descModel] = getSpacialHistogramDescriptors(ptsModel, sample_ptsModel, min_pts, R, thVar, ALIGN_POINTS, CENTER, k);
    if RAND_CROP
        [featRand, descRand] = getSpacialHistogramDescriptors(ptsRand, sample_ptsRand, min_pts, R, thVar, ALIGN_POINTS, CENTER, k);
    end
    [featSurface, descSurface] = getSpacialHistogramDescriptors(ptsSurface, sample_ptsSurface, min_pts, R, thVar, ALIGN_POINTS, CENTER, k);        
else
    error('Descriptor method must be either Moment or Rotational');
end

    
%% Apply weights to descriptors
if strcmp(descriptor, 'Moment')    
    weights.M1 = 0.3; % 0
    weights.M2 = 0.75; % 0.3
    weights.M3 = 0.38; % 0.8
    weights.M4 = 0.54; % 0.5

    descSurfaceW = applyMomentWeight(descSurface, weights);
    descModelW = applyMomentWeight(descModel, weights);
    if RAND_CROP
        descRandW = applyMomentWeight(descRand, weights);
    end
else
    descSurfaceW = descSurface;
    descModelW = descModel;
    if RAND_CROP
        descRandW = descRand;
    end
end
    
%% Modify descriptors based on their variance --> Mahalanobis Distance
% get covariance matrix of all descriptors (this is a 31x31 matrix)
S_Model = cov([descModelW; descSurfaceW]);
S_Rand = cov([descRandW; descSurfaceW]);

% get square root, so that A'*A = S (S is hermitian)
A_Model = real(S_Model^0.5);
A_Rand = real(S_Rand^0.5);

% Modify Descriptors accordingly
descModel_Mnobis = descModelW*A_Model;
descSurfaceM_Mnobis = descSurfaceW*A_Model;

descRand_Mnobis = descRandW*A_Rand;
descSurfaceR_Mnobis = descSurfaceW*A_Rand;

%% Match features between Surface and Model / Random Crop

% Define matching algorithm parameters
par.MAHALANOBIS = false; % use Mahalanobis Distance instead of L2
par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold =  0.5; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.6; % 0.6 (default) nearest neighbor ambiguity rejection
par.Metric =  'SSD'; % SSD (default) for L2, SAD for L1
par.Unique = false; % true: 1-to-1 mapping only, else set false (default)

if par.MAHALANOBIS
    dSM = descSurfaceM_Mnobis;
    dMM = descModel_Mnobis;
    dSR = descSurfaceR_Mnobis;
    dRR = descRand_Mnobis;
else
    dSM = descSurfaceW;
    dMM = descModelW;
    dSR = descSurfaceW;
    dRR = descRandW;
end

matchesModel = matchFeatures(dSM, dMM, ...
        'Method', par.Method, ...
        'MatchThreshold', par.MatchThreshold, ... 
        'MaxRatio', par.MaxRatio, ... 
        'Metric', par.Metric, ...
        'Unique', par.Unique); 
if RAND_CROP
    matchesRand = matchFeatures(dSR, dRR, ...
            'Method', par.Method, ...
            'MatchThreshold', par.MatchThreshold, ... 
            'MaxRatio', par.MaxRatio, ... 
            'Metric', par.Metric, ...
            'Unique', par.Unique); 
end
    
%% Get distance between matching points
% this makes sense, because the pointclouds are already aligned. The
% distance between matching INLIERS of the model will thus be small and a
% simple distance threshold can be used to determine whether a match is an
% inlier 

% maxDist specifies matching distance to count inliers 
maxDist = 1.3; 

% matches of surface and model
loc1S = featSurface(matchesModel(:, 1), :);
loc1M = featModel(matchesModel(:, 2), :);
d1 = vecnorm(loc1M - loc1S, 2, 2);
inlier_idx = find(d1 < maxDist);
inliers1 = length(inlier_idx);

% precision
inliersPrecision = inliers1/size(d1, 1)*100;

if RAND_CROP
    % matches of surface and random crop
    loc2S = featSurface(matchesRand(:, 1), :);
    loc2M = featRand(matchesRand(:, 2), :);
    d2 = vecnorm(loc2M - loc2S, 2, 2);
    inliers2 = length(find(d2 <= maxDist));
end

%% helper function that applies weight to descriptor
function descW = applyMomentWeight(desc, weights)
    descW = desc;
    descW(:, 1:6) = desc(:, 1:6)*weights.M1;
    descW(:, 7:12) = desc(:, 7:12)*weights.M2;
    descW(:, 13:22) = desc(:, 13:22)*weights.M3;
    descW(:, 23:31) = desc(:, 23:31)*weights.M4;
end

%% helper function: center PC (0,0) in middle of crop
function pc_out = centerPointCloud(pc)
    % translation of pointcloud
    T = [(pc.XLimits(1) + pc.XLimits(2))/2, ...
        (pc.YLimits(1) + pc.YLimits(2))/2, ...
        (pc.ZLimits(1) + pc.ZLimits(2))/2];
    A = eye(4);
    A(4, 1:3) = -T;
    tform = affine3d(A);
    pc_out = pctransform(pc, tform);
end

%% helper function: uniformly sample point cloud
function sample_pts = pcUniformSamples(pcIn, d) 
    % create meshgrid of points	  
    x = pcIn.XLimits(1):d:pcIn.XLimits(2);
    y = pcIn.YLimits(1):d:pcIn.YLimits(2);
    z = pcIn.ZLimits(1):d:pcIn.ZLimits(2);
    [X,Y,Z] = meshgrid(x,y,z);
    sample_pts = cat(4, X, Y, Z);
    sample_pts = reshape(sample_pts, [], 3);
end

%% helper function: random-uniformly sample point cloud
function sample_pts = pcRandomUniformSamples(pcIn, d)
    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcIn.XLimits(2) - pcIn.XLimits(1);
    rangeY = pcIn.YLimits(2) - pcIn.YLimits(1);
    rangeZ = pcIn.ZLimits(2) - pcIn.ZLimits(1);
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers
    sample_pts = rand(num_pts, 3);
    
    % scale numbers so that they fit into the correct range
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    sample_pts = sample_pts + ...
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)];
end