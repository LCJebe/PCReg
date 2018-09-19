%% completeExperimentFast.m
% running complete experiment, in a sparse and fast version
rng('shuffle');

%% options

%% setup
centerSurface = [34.8970   21.9820   25.4546];

%% read Model and Surface pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothUp3.pcd'));
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));
pcSurface = centerPointCloud(pcSurface); % this has been used for the descriptors

%% read descriptors for model and surface
% _raw --> SurfaceNew_DS3.pcd
% _rand --> SurfaceNew_DS3.pcd transformed by TF (saved variable)
% _rand2 --> SurfaceNew_DS3_rand2.pcd

load('Data/Descriptors/featSurface0.3_raw.mat');
load('Data/Descriptors/descSurface0.3_raw.mat'); 
load('Data/Descriptors/featModel0.3.mat');
load('Data/Descriptors/descModel0.3.mat'); 

    
%% sliding sphere selection of descriptor locations
% sparse, because this is the FAST experiment
% setup: get length of cuboid diagonal
pcSurfaceRaw = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));
xLimS = pcSurfaceRaw.XLimits;
yLimS = pcSurfaceRaw.YLimits;
zLimS = pcSurfaceRaw.ZLimits;
diagSurface = norm([xLimS(2) - xLimS(1), ...
                    yLimS(2) - yLimS(1), ...
                    zLimS(2) - zLimS(1)]);
                
% "smaller" setting, even though we have -3.5, because descriptors are
% selected within radius of R_desc now, and before we cropped with R_crop
% and selected descriptors within R_crop - 3.5.
% so: R_desc = R_crop - 3.5.
R_desc = diagSurface / 2 - 3.5; 


% make pointcloud out of descriptor locations
pcDescModel = pointCloud(featModel);

% sample uniform grid again for large spheres in which we are matching
d_Spheres = 5;
sample_ptsSpheres = pcUniformSamples(pcDescModel, d_Spheres);

% first: check for valid spheres by number of descriptors in sphere
tic
min_pts = 1400; % high threshold, to be faster
max_pts = inf;
valid_mask = false(size(sample_ptsSpheres, 1), 1);
num_desc = zeros(size(sample_ptsSpheres, 1), 1);
parfor i = 1:size(sample_ptsSpheres, 1)
    c = sample_ptsSpheres(i, :)
    [~, dists] = getLocalPoints(featModel, R_desc, c, min_pts, max_pts);
    if ~isempty(dists)
        valid_mask(i) = 1;
        num_desc(i) = length(dists);
    end
end
fprintf('Obtained %d test positions in %0.1f seconds...\n', sum(valid_mask), toc);

% retain only valid sample_pts
num_desc = num_desc(valid_mask);
valid_mask = cat(2, valid_mask, valid_mask, valid_mask);
sample_ptsSpheres = reshape(sample_ptsSpheres(valid_mask), [], 3);

%% match the points and get number of putative matches for each sphere
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

par.VERBOSE = 1;

num_spheres = size(sample_ptsSpheres, 1);
num_putative_all = zeros(num_spheres, 1);
num_desc_all = zeros(num_spheres, 1);
matchesCells = cell(num_spheres, 1);
featureCells = cell(num_spheres, 1);
locationArray = nan(num_spheres, 3);

tic
fprintf('Getting Putative Matches for all spheres...\n');
fprintf('This will take about %0.1f x t minutes...\n', double(num_spheres)/4/60*1.5);

% do some calculations before the parfor to avoid running out of memory
portion_size = 16;
num_portions = ceil(double(num_spheres)  / portion_size);
for portion = 1:num_portions

    i_start = 1+(portion-1)*portion_size;
    
    featCur_portion = cell(portion_size, 1);
    descCur_portion = cell(portion_size, 1);
    locationCur_portion = nan(portion_size, 3);
    
    for i = i_start:i_start+portion_size-1 % counter through all spheres
        if i > num_spheres
            break;
        end
        
        ii = i - i_start + 1; % counter from 1 to portion_size
        c = sample_ptsSpheres(i, :);

        % mask of descriptors that lie within sphere
        mask = getDescriptorMask(featModel, c, R_desc, 0);

        % based on the mask, return the relevant sample points and descriptors
        maskF = cat(2, mask, mask, mask);
        featCur = reshape(featModel(maskF), [], 3);  
        maskD = repmat(mask, 1, size(descModel, 2));
        descCur = reshape(descModel(maskD), [], size(descModel, 2));
        
        featCur_portion{ii} = featCur;
        descCur_portion{ii} = descCur;
        locationCur_portion(ii, :) = c;
    end

    endindex = min(i_start+portion_size-1, num_spheres);
    parfor i = i_start:endindex
        ii = i - i_start + 1;
        descCur = descCur_portion{ii};
        featCur = featCur_portion{ii};
        c = locationCur_portion(ii, :);
        
        % now match with the specified region
        matchesModel = getMatches(descSurface, descCur, par);

        % save number of putative matches
        num_putative_all(i) = size(matchesModel, 1);
        num_desc_all(i) = size(descCur, 1);
        matchesCells{i} = matchesModel; % indices of matches
        featureCells{i} = featCur; % locations of local features
        locationArray(i, :) = c; % center of large sphere
    end
end

%% reduce the cell arrays to contain only non-empty cells
empty_idx = cellfun('isempty',featureCells);
matchesCellsFull = matchesCells(~empty_idx);
featureCellsFull = featureCells(~empty_idx);
locationArrayFull = locationArray(~empty_idx, :);
num_putativeFull = num_putative_all(~empty_idx);
num_descFull = num_desc_all(~empty_idx);

fprintf('Matched to all positions in %0.1f seconds...\n', toc);

% clear some memory
clear featureCells matchesCells locationArray 

%% perform RANSAC on spheres with a lot of putative matches
putative_thresh = 170;

% RANSAC options
options.minPtNum = 3; % 3
options.iterNum = 1e4; % 3e4
options.thDist = 0.3;  % 0.2
options.thInlrRatio = 0.08; % 0.1
options.REFINE = true;
options.VERBOSE = 1;

% apply threshold
idx = find(num_putativeFull>putative_thresh);

% prepare cells to contain trial cases only
matchesCellsTrial = matchesCellsFull(idx);
featureCellsTrial = featureCellsFull(idx);
locationArrayTrial = locationArrayFull(idx, :);
num_putativeTrial = num_putativeFull(idx);
num_descTrial = num_descFull(idx);


num_trials = length(idx);


% little print for timing
fprintf('Performing RANSAC on all %d promising spheres...\n', num_trials);
fprintf('This will take about %0.1f x t minutes...\n', double(num_trials)/4/60);

% statistics: numPutative, numSuccess, maxInliers, maxInlierRatio
statsPutative = zeros(num_trials, 1);
statsSuccess = zeros(num_trials, 1);
statsInliers = zeros(num_trials, 1);
statsRatio = zeros(num_trials, 1);
transformCells = cell(num_trials, 1);

parfor i = 1:num_trials
    
    matches = matchesCellsTrial{i};
    featuresM = featureCellsTrial{i};
    pts1 = featSurface(matches(:, 1), :);
    pts2 = featuresM(matches(:, 2), :);
    
    % sanity check
    assert(num_putativeTrial(i) == size(matches, 1));
    
    % RANSCAC
    tic
    [T, inlierPtIdx, numSuccess, maxInliers, maxInlierRatio] = ...
                ransac(pts1,pts2,options,@estimateTransform,@calcDists);
    toc
           
    % save statistics. needed to find out which sphere matches best
    
    statsNumDesc(i) = num_descTrial(i);
    statsPutative(i) = num_putativeTrial(i);
    statsSuccess(i) = numSuccess;
    statsInliers(i) = maxInliers;
    statsRatio(i) = maxInlierRatio;
    transformCells{i} = T;
end

%%% EVALUTATION %%%

%% select points that fulfill certrain criteria
c = centerSurface;

%threshold for green points
thSucc = 0;
thInliers = 28;
thRatio = 10; % in percent
thPutative = 170;

mask = statsSuccess >= thSucc & statsInliers >= thInliers & statsRatio >= thRatio & statsPutative >= thPutative;
indices = find(mask > 0);

% get location and stats of only successful runs
statsPutativeGood = statsPutative(indices);
%statsNumDescGood = statsNumDesc(indices)';
statsSuccessGood = statsSuccess(indices);
statsInliersGood = statsInliers(indices);
statsRatioGood = statsRatio(indices);
goodLocations = locationArrayTrial(indices, :);
transformCellsGood = transformCells(indices);

% show all points in green and the correct center in red
goodLocations = [goodLocations; centerSurface];
color = ones(size(goodLocations, 1), 3) .* [0, 255, 0];
color(end, :) = [255, 0, 0];

% show all sampled sphere positions in black
SHOW_SAMPLES = true;
if SHOW_SAMPLES
    good_samples = getLocalPoints(locationArrayTrial, inf, centerSurface, 0, inf)+centerSurface;
    showLocations = [good_samples; goodLocations];
    color =  [ones(size(good_samples, 1), 3).*[0.4, 0.4, 0.4]; color];
end

pcshow(pointCloud(showLocations, 'Color', color), 'MarkerSize', 50);

%% finally, cluster points
r = 1.6*d_Spheres; % 1.6 = sqrt(3) because of diagonal adjacent points
clusters = clusterPoints(goodLocations(1:end-1, :), r);
cluster_size = cellfun('length', clusters);

%% run RANSAC once for each cluster
% with GLOBAL Alignment! That is, the surface descriptors have to be
% calculated again according to the previously found RANSAC rotation

num_clusters = length(clusters);
finalMatchCells = cell(num_clusters, 1);
finalFeatModelCells = cell(num_clusters, 1);
finalFeatSurfaceCells = cell(num_clusters, 1);
precision_test = zeros(num_clusters, 1);

for i = 1:num_clusters    
    % choose point in the middle of cluster
    clusterPointsCur = goodLocations(clusters{i}, :);
    locCur = mean(clusterPointsCur, 1);
    
    % get transform from that cluster (closest point to mean)
    dists = vecnorm(clusterPointsCur - locCur, 2, 2);
    [~, idx] = min(dists);
    transCur = transformCellsGood{clusters{i}(idx)};
    
    % transform the Surface accordingly (like autoalign)
    pts_Surface_tform = quickTF(pcSurface.Location, invertTF(transCur));
    
    % get new descriptors from the surface, without local alignment
    d = 0.5; 
    margin = 3.5;
    sample_ptsSurface = pcRandomUniformSamples(pointCloud(pts_Surface_tform), d, margin);

    % descriptor options
    descOpt.ALIGN_POINTS = false; % this false is the key
    descOpt.min_pts = 500;
    descOpt.max_pts = 6000;
    descOpt.R = 3.5; % should equal margin
    descOpt.thVar = [3, 1.5]; 
    descOpt.k = 0.85;
    descOpt.VERBOSE = 1;
    descOpt.max_region_size = 15;

    [featSurface_noLRF, descSurface_noLRF] = ...
            getSpacialHistogramDescriptors(pts_Surface_tform, sample_ptsSurface, descOpt);
        
    % get features and descriptors from the Model as well (no LRF)
    featModel_noLRF = load('Data/Descriptors/featModel0.3_noLRF.mat');
    descModel_noLRF = load('Data/Descriptors/descModel0.3_noLRF.mat');
    featModel_noLRF = featModel_noLRF.featModel;
    descModel_noLRF = descModel_noLRF.descModel;
    
    % select correct spherical subest of descriptors based on locCur
    mask = getDescriptorMask(featModel_noLRF, locCur, R_desc, 0);
    indices = find(mask);
    featCur_noLRF = featModel_noLRF(mask, :);
    descCur_noLRF = descModel_noLRF(mask, :);
    
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

    par.VERBOSE = 1;

    % match
    tic
    matchesCur = getMatches(descSurface_noLRF, descCur_noLRF, par);
    
    %save matches
    finalMatchCells{i} = matchesCur;
    finalFeatModelCells{i} = featCur_noLRF;
    finalFeatSurfaceCells{i} = featSurface_noLRF;
end

%% get inliers for each tested cluster by distance check

inlier_precisions = zeros(num_clusters, 1);
for i = 1:num_clusters
    matchesCur = finalMatchCells{i};
    featCur_noLRF = finalFeatModelCells{i};
    featSurface_noLRF = finalFeatSurfaceCells{i};

    % Perform distance check on matches    
    maxDist = 1.5; 

    % matches of surface and model
    pts1 = featSurface_noLRF(matchesCur(:, 1), :);
    pts2 = featCur_noLRF(matchesCur(:, 2), :);
    d1 = vecnorm(pts1 - pts2, 2, 2);
    inlier_idx = find(d1 < maxDist);
    inliers1 = length(inlier_idx);

    % precision
    inliersPrecision = inliers1/size(d1, 1)*100;    
    
    % add inlier ratio to later find out the correct cluster
    inlier_precisions(i) = inliersPrecision;
end

%% use the best hit with the most inliers for a refined transformation estimation
[maxInliers, bestCluster] = max(inlier_precisions);

% get inliers again
matchesCur = finalMatchCells{bestCluster};
featCur_noLRF = finalFeatModelCells{bestCluster};
featSurface_noLRF = finalFeatSurfaceCells{bestCluster};

% Perform distance check on matches    
maxDist = 1.5; 

% matches of surface and model
pts1 = featSurface_noLRF(matchesCur(:, 1), :);
pts2 = featCur_noLRF(matchesCur(:, 2), :);
d1 = vecnorm(pts1 - pts2, 2, 2);
inlier_idx = find(d1 < maxDist);

% should be close tot eye(4) because it's a refinement
T_refine = estimateTransform(pts1(inlier_idx, :), pts2(inlier_idx, :));

% apply refinement to surface, and save as final surface
pts_Surface_final = quickTF(pts_Surface_tform, invertTF(T_refine));

% create colorful pointcloud
color = pcSurface.Color;
pcSurface_final = pointCloud(pts_Surface_final, 'Color', color);

% save
path = 'Data/PointClouds/';
save_file = strcat(path, 'Surface_final.pcd');
pcwrite(pcSurface_aligned, save_file);


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

%% helper function: descriptors within certain distance of given center point
function mask = getDescriptorMask(featModel, current_center, R_desc, margin)
    feat_rel = featModel - current_center;
    dists = vecnorm(feat_rel, 2, 2);
    mask = dists < (R_desc + margin); % note that margin < 0
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