%% slideMatchingWindow_v2.m
% experiment, where we slowly slide the crop from the ideal position
% outwards with a ramdom angle, and record how well RANSAC performs with
% different metrics (inlier percentage, num successes, num inliers)
% DIFFERENCE TO v1:
% this script uses pre-caluculated descriptors on the whole model

%% load limits of aligned surface for crop-reference!
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3_alignedM.pcd'));
xLimS = pcSurface.XLimits;
yLimS = pcSurface.YLimits;
zLimS = pcSurface.ZLimits;

% get center of aligned Surface
centerSurface = [(xLimS(2) + xLimS(1))/2, ...
                 (yLimS(2) + yLimS(1))/2, ...
                 (zLimS(2) + zLimS(1))/2];
             
% for cuboid diagonal, use Surface that is aligned with the axes, not model
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));
xLimS = pcSurface.XLimits;
yLimS = pcSurface.YLimits;
zLimS = pcSurface.ZLimits;

% get length of cuboid diagonal
diagSurface = norm([xLimS(2) - xLimS(1), ...
                    yLimS(2) - yLimS(1), ...
                    zLimS(2) - zLimS(1)]);
                
%% transform Surface randomly and center
% transform surface manually to random position and orientation
MANUAL_TF = true;
RotS = rand(1, 3)*2*pi;
transS = [13, 25, -17];
pcSurface = pcRigidBodyTF(pcSurface, RotS, transS);

% shift pointcloud to center
pcSurface = centerPointCloud(pcSurface);

%% load pre-calculated descriptors
load('featModel0.4.mat');
load('descModel0.4.mat');


%% specify crop area, used to crop out descriptors by location
R_crop = diagSurface/2; % "smaller"

% now define a random direction of sliding
% optional: load experiment
LOAD_EXPERIMENT = false;
if LOAD_EXPERIMENT
    load('putative_curve_mean.mat');
    load('success_curve_mean.mat');
    load('inlier_curve_mean.mat');
    load('ratio_curve_mean.mat');
    load('ex.mat');
    start_ex = ex;
end
num_experiments = 10;
for ex = start_ex:num_experiments+start_ex - 1
    rng('shuffle');
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
    P = pcModel.Location; % (N x 3)
    SE = E-S; % (1 x 3);
    PS = P-S; % (N x 3)
    dQS = (PS*SE')/(SE*SE'); % (N x 1)
    QS =  dQS.* SE; % (N x 3)
    d = vecnorm(QS, 2, 2); % (N x 1)
    mask = (d < norm(SE)) & (vecnorm(QS-PS, 2, 2) < R_crop) & (QS*SE' > 0); % mask for points inside crop region
    mask = cat(2, mask, mask, mask);
    pts_crop = reshape(P(mask), [], 3);  

    color_crop = reshape(pcModel.Color(mask), [], 3);  

    % show cropped region
    % pcshow(pointCloud(pts_crop, 'Color', color_crop));

    %% now calculate descriotors on the cropped region 
    pcModel = pointCloud(pts_crop);

    % define steps between 0 and 1 which let the sphere slide from
    % centerSurface (0, 0, 0) to centerOff = 2*R_crop*[x, y, z]

    steps = 0:0.1:1;

    % get sampling points for the model and the crop
    dS = 0.3;
    dM = 0.4;
    margin = 3.5;
    sample_ptsSurface = pcRandomUniformSamples(pcSurface, dS, margin);
    sample_ptsModel = pcCylindricalSamples(pcModel, dM, margin, S, E, R_crop);

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
    [featModel, descModel] = ...
            getSpacialHistogramDescriptors(pcModel.Location, sample_ptsModel, descOptM);

    %% Now run the experiment 
    % define steps between 0 and 1 which let the sphere slide from
    % centerSurface  to centerOff 
    steps = 0:0.05:1;

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
    options.iterNum = 2e4; 
    options.thDist = 0.5; 
    options.thInlrRatio = 0.1; 
    options.REFINE = true;

    % minimum number of putative matches
    min_matches = 50;

    % placeholder for the plots
    putative_curve = zeros(length(steps), 1); 
    success_curve = zeros(length(steps), 1);
    inlier_curve = zeros(length(steps), 1);
    ratio_curve = zeros(length(steps), 1);

    for i = 1:length(steps)
        step = steps(i);
        current_center = centerSurface + step*(centerOff - centerSurface);
        mask = getDescriptorIndices(featModel, current_center, R_crop, -3.5);

        % based on the mask, return the relevant sample points and descriptors
        maskF = cat(2, mask, mask, mask);
        featCur = reshape(featModel(maskF), [], 3);  
        maskD = repmat(mask, 1, size(descModel, 2));
        descCur = reshape(descModel(maskD), [], size(descModel, 2));

        % now match with the specified region
        matchesModel = getMatches(descSurface, descCur, par);

        % run RANSAC and return a few metrics
        pts1 = featSurface(matchesModel(:, 1), :);
        pts2 = featCur(matchesModel(:, 2), :);
        if size(pts1, 1) > min_matches
            [T, inlierPtIdx, numSuccess, maxInliers, maxInlierRatio] = ...
                ransac(pts1,pts2,options,@estimateTransform,@calcDists, options.REFINE);
        else
            [T, inlierPtIdx, numSuccess, maxInliers, maxInlierRatio] = deal(0);
        end

        putative_curve(i) = size(matchesModel, 1);
        success_curve(i) = numSuccess;
        inlier_curve(i) = maxInliers;
        ratio_curve(i) = maxInlierRatio;
    end
    
    % iterative mean of curves
    if ex == 1
        putative_curve_mean = putative_curve;
        success_curve_mean = success_curve;
        inlier_curve_mean = inlier_curve;
        ratio_curve_mean = ratio_curve;
    else
        putative_curve_mean = putative_curve_mean + 1/(ex+1)*(putative_curve-putative_curve_mean);
        success_curve_mean = success_curve + 1/(ex+1)*(success_curve-success_curve_mean);
        inlier_curve_mean = inlier_curve + 1/(ex+1)*(inlier_curve-inlier_curve_mean);
        ratio_curve_mean = ratio_curve + 1/(ex+1)*(ratio_curve-ratio_curve_mean);
    end
    
    % save to workspace
    if ~mod(ex, 10)
        save('putative_curve_mean.mat', 'putative_curve_mean');
        save('success_curve_mean.mat', 'success_curve_mean');
        save('inlier_curve_mean.mat', 'inlier_curve_mean');
        save('ratio_curve_mean.mat', 'ratio_curve_mean');
        save('ex.mat', 'ex');
    end
    
end % for ex
% save to workspace
save('putative_curve_mean.mat', 'putative_curve_mean');
save('success_curve_mean.mat', 'success_curve_mean');
save('inlier_curve_mean.mat', 'inlier_curve_mean');
save('ratio_curve_mean.mat', 'ratio_curve_mean');
save('ex.mat', 'ex');


%% Plots
close all;
figtitle = sprintf('Theta: %0.5f, Phi: %0.5f', theta, phi);
figure('Name', figtitle);
subplot(4, 1, 1);
plot(steps, putative_curve);
title('Putative Matches');
grid;
subplot(4, 1, 2);
plot(steps, success_curve);
title('Number of Successes');
grid;
subplot(4, 1, 3);
plot(steps, inlier_curve);
title('Number of Inliers');
grid;
subplot(4, 1, 4);
plot(steps, ratio_curve);
title('Ratio of Inliers');
grid;

%% helper function: descriptors within certain distance of given center point
function mask = getDescriptorIndices(featModel, current_center, R_crop, margin)
    feat_rel = featModel - current_center;
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
        [pcIn.XLimits(1), pcIn.YLimits(1), pcIn.ZLimits(1)] - margin;
    
    % sample points are called "P" for shortness
    % SE: vector from Startpoint to Endpoint of Cylinder
    SE = E-S; % (1 x 3);
    % Vector from each point to the starting point
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