%% slideMatchingWindow.m
% experiment, where we slowly slide the crop from the ideal position
% outwards with a ramdom angle, and record how well RANSAC performs with
% different metrics (inlier percentage, num successes, num inliers)

%% read Model pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothUp3.pcd'));
pcSurfce = pcread(strcat(path, 'Surface_DS3.pcd'));

% load limits of aligned surface for crop-reference!
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3_alignedM.pcd'));
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


%% specify crop area
R_crop = diagSurface/2; % smaller

% now define a random direction of sliding
theta = 2*pi*rand(1, 1);
phi = acos(2*rand(1, 1)-1);

% convert to x, y, z
x = sin(theta)*cos(phi); % 0 for theta = 0, phi = 0
y = sin(theta)*sin(phi); % 0 for theta = 0, phi = 0
z = cos(theta);          % 1 for theta = 0, phi = 0

% complete sliding offset will be 2*R_crop
centerOff = centerSurface + 2*R_crop*[x, y, z];

% get start and endpoints
S = centerSurface - R_crop*[x, y, z];
E = centerSurface + 3*R_crop*[x, y, z];

% example for list of points P of size (N x 3)
%P = pcModel.Location; % (N x 3)
P = pcModel.Location;
SE = E-S; % (1 x 3);
PS = P-S; % (N x 3)
dQS = (PS*SE')/(SE*SE'); % (N x 1)
QS =  dQS.* SE; % (N x 3)
d = vecnorm(QS, 2, 2); % (N x 1)
mask = (d < norm(SE)) & (vecnorm(QS-PS, 2, 2) < R_crop) & (QS*SE' > 0); % mask for points inside crop region
mask = cat(2, mask, mask, mask);
pts_crop = reshape(P(mask), [], 3);  

% show cropped region
%pcshow(pts_crop);

%% now calculate descriotors on the cropped region 
pcCrop = pointCloud(pts_crop);

% define steps between 0 and 1 which let the sphere slide from
% centerSurface (0, 0, 0) to centerOff = 2*R_crop*[x, y, z]

steps = 0:0.1:1;

% get sampling points for the model and the crop
d = 0.4;
margin = 3.5;
sample_ptsSurface = pcRandomUniformSamples(pcSurface, d, margin);
sample_ptsCrop = pcCylindricalSamples(pcCrop, d, margin, S, E, R_crop);

% descriptor options
descOptM.ALIGN_POINTS = true;
descOptM.CENTER = false;
descOptM.min_pts = 500;
descOptM.max_pts = 6000;
descOptM.R = 3.5;
descOptM.thVar = [3, 1.5]; 
descOptM.k = 0.85;

% get descriptors and locations
[featSurface, descSurface] = ...
        getSpacialHistogramDescriptors(pcSurface.Location, sample_ptsSurface, descOptM);   
[featCrop, descCrop] = ...
        getSpacialHistogramDescriptors(pcCrop.Location, sample_ptsCrop, descOptM);
    
%% Now run the experiment 
% define steps between 0 and 1 which let the sphere slide from
% centerSurface  to centerOff 
steps = 0:0.1:1;

% matching options
par.UNNORMALIZE = true;
par.norm_factor = 2; % 2

par.CHANGE_METRIC = true;
par.metric_factor = 0.6; % 0.45 / 0.6

par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold = 10; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.99; % 0.6 (default) nearest neighbor ambiguity rejection
par.Metric =  'SAD'; % SSD (default) for L2, SAD for L1
par.Unique = true; % true: 1-to-1 mapping only, else set false (default)

% RANSAC options
options.minPtNum = 3; 
options.iterNum = 5e3; 
options.thDist = 0.5; 
options.thInlrRatio = 0.1; 
options.REFINE = true;

% placeholder for the plots
success_curve = zeros(length(steps), 1);
inlier_curve = zeros(length(steps), 1);
ratio_curve = zeros(length(steps), 1);

for i = 1:length(steps)
    step = steps(i);
    current_center = centerSurface + step*(centerOff - centerSurface);
    mask = getDescriptorIndices(featCrop, current_center, R_crop, -3.5);
        
    % based on the mask, return the relevant sample points and descriptors
    maskF = cat(2, mask, mask, mask);
    featCur = reshape(featCrop(maskF), [], 3);  
    maskD = repmat(mask, 1, size(descCrop, 2));
    descCur = reshape(descCrop(maskD), [], size(descCrop, 2));
    
    % now match with the specified region
    matchesCur = getMatches(descSurface, descCur, par);
    
    % run RANSAC and return a few metrics
    pts1 = featSurface;
    pts2 = featCur;
    [T, inlierPtIdx, numSuccess, maxInliers, maxInlierRatio] = ...
        ransac(pts1,pts2,options,@estimateTransform,@calcDists, options.REFINE);
    
    success_curve(i) = numSuccess;
    inlier_curve(i) = maxInliers;
    ratio_curve(i) = maxInlierRatio;
end

%% Plots
close all;
figure(1);
subplot(3, 1, 1);
plot(steps, success_curve);
title('Number of Successes');
subplot(3, 1, 2);
plot(steps, inlier_curve);
title('Number of Inliers');
subplot(3, 1, 3);
plot(steps, ratio_curve);
title('Ratio of Inliers');

%% helper function: descriptors within certain distance of given center point
function mask = getDescriptorIndices(featCrop, current_center, R_crop, margin)
    feat_rel = featCrop - current_center;
    dists = vecnorm(feat_rel, 2, 2);
    mask = dists < (R_crop + margin); % note that margin < 0
end


%% helper function: sample points for spherical crop
function sample_pts = pcCylindricalSamples(pcIn, d, margin, S, E, R_crop)

    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcIn.XLimits(2) - pcIn.XLimits(1) + 2*margin;
    rangeY = pcIn.YLimits(2) - pcIn.YLimits(1) + 2*margin;
    rangeZ = pcIn.ZLimits(2) - pcIn.ZLimits(1) + 2*margin;
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers in range [0, 1]
    sample_pts = rand(num_pts, 3); 
    
    % scale numbers so that they fit into the correct range, without using
    % margin yet
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    P = sample_pts + ...
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)];
    
    % sample points are called "P" for shortness
    SE = E-S; % (1 x 3);
    PS = P-S; % (N x 3)
    dQS = (PS*SE')/(SE*SE'); % (N x 1)
    QS =  dQS.* SE; % (N x 3)
    d = vecnorm(QS, 2, 2); % (N x 1)
    mask = (d < norm(SE)) & (vecnorm(QS-PS, 2, 2) < (R_crop-margin)) & (QS*SE' > 0); % mask for points inside crop region
    mask = cat(2, mask, mask, mask);
    sample_pts = reshape(P(mask), [], 3);  
end

%% helper function: random-uniformly sample point cloud
function sample_pts = pcRandomUniformSamples(pcIn, d, margin)
    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcIn.XLimits(2) - pcIn.XLimits(1) + 2*margin;
    rangeY = pcIn.YLimits(2) - pcIn.YLimits(1) + 2*margin;
    rangeZ = pcIn.ZLimits(2) - pcIn.ZLimits(1) + 2*margin;
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers in range [0, 1]
    sample_pts = rand(num_pts, 3); 
    
    % scale numbers so that they fit into the correct range
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    sample_pts = sample_pts + ...
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)] - margin;
end


%% function that calculates the distance between points after transform T
% used as distance metric for RANSAC
function d = calcDists(T,pts1,pts2)
    %	Project PTS2 to PTS2_trans using the rigid transform T, then calcultate the distances between
    %	PTS1 and PTS2_trans

    pts2(:, 4) = 1;
    pts2_trans = pts2*T;
    pts2_trans = pts2_trans(:, 1:3);
    d = sum((pts1-pts2_trans).^2,2);
end