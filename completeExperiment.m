%% slideMatchingWindow.m
% experiment, where we slowly slide the crop from the ideal position
% outwards with a ramdom angle, and record how well RANSAC performs with
% different metrics (inlier percentage, num successes, num inliers)
rng('shuffle');

%% options
% load descriptors from folder instead of calculating them again
LOAD_DATA = true;


%% read Model and Surface pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothUp3.pcd'));
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));

if ~LOAD_DATA
    %% get descriptors for surface and for the complete model
    dS = 0.3; % denser (0.2 or 0.3)
    sampleOptM.d = 0.4; % sparser (0.4 or more)
    sampleOptM.margin = 3.5; % should equal R

    sample_ptsSurface = pcRandomUniformSamples(pcSurface, dS, sampleOptM.margin);

    % descriptor options
    descOptM.ALIGN_POINTS = true;
    descOptM.CENTER = false;
    descOptM.min_pts = 500;
    descOptM.max_pts = 6000;
    descOptM.R = 3.5; % should equal margin
    descOptM.thVar = [3, 1.5]; 
    descOptM.k = 0.85;
    descOptM.VERBOSE = 1;
    descOptM.max_region_size = 15;

    descOptS = descOptM;
    descOptS.max_region_size = 100;


    [featModel, descModel] = ...
            speedyDescriptors(pcModel.Location, sampleOptM, descOptM);   
    % save the features and descriptors to workspace
    save('Data/Descriptors/featModel0.4.mat', 'featModel');
    save('Data/Descriptors/descModel0.4.mat', 'descModel');

    [featSurface, descSurface] = ...
            getSpacialHistogramDescriptors(pcSurface.Location, sample_ptsSurface, descOptS);   
else
    load('Data/Descriptors/featSurface0.3.mat');
    load('Data/Descriptors/descSurface0.3.mat'); 
    load('Data/Descriptors/featModel0.3.mat');
    load('Data/Descriptors/descModel0.3.mat'); 
end
    
%% sliding sphere selection of descriptor locations
% setup: get length of cuboid diagonal
xLimS = pcSurface.XLimits;
yLimS = pcSurface.YLimits;
zLimS = pcSurface.ZLimits;
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
d_Spheres = 2;
sample_ptsSpheres = pcUniformSamples(pcDescModel, d_Spheres);

% first: check for valid spheres by number of descriptors in sphere
tic
min_pts = 200;
max_pts = inf;
valid_mask = false(size(sample_ptsSpheres, 1), 1);
parfor i = 1:size(sample_ptsSpheres, 1)
    c = sample_ptsSpheres(i, :);
    [~, dists] = getLocalPoints(featModel, R_desc, c, min_pts, max_pts);
    if ~isempty(dists)
        valid_mask(i) = 1;
    end
end
fprintf('Obtained %d test positions in %0.1f seconds...\n', sum(valid_mask), toc);

% retain only valid sample_pts
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

num_putative_all = zeros(size(sample_ptsSpheres, 1), 1);
matchesCells = cell(size(sample_ptsSpheres, 1), 1);
featureCells = cell(size(sample_ptsSpheres, 1), 1);
locationArray = nan(size(sample_ptsSpheres, 1), 3);

tic
num_spheres = size(sample_ptsSpheres, 1);
fprintf('Getting Putative Matches for all spheres...\n');
fprintf('This will take about %d x t minutes...\n', round(num_spheres/4/60));

% do some calculations before the parfor to avoid running out of memory
portion_size = 16;
num_portions = ceil(double(num_spheres)  / portion_size);
for portion = 1:num_portions

    i_start = 1+(portion-1)*portion_size;
    
    featCur_portion = cell(portion_size, 1);
    descCur_portion = cell(portion_size, 1);
    
    for i = i_start:i_start+portion_size-1 % counter through all spheres
        ii = i - i_start + 1; % counter from 1 to portion_size
        c = sample_ptsSpheres(i, :);

        % mask of descriptors that lie within sphere
        mask = getDescriptorIndices(featModel, c, R_desc, 0);

        % based on the mask, return the relevant sample points and descriptors
        maskF = cat(2, mask, mask, mask);
        featCur = reshape(featModel(maskF), [], 3);  
        maskD = repmat(mask, 1, size(descModel, 2));
        descCur = reshape(descModel(maskD), [], size(descModel, 2));
        
        featCur_portion{ii} = featCur;
        descCur_portion{ii} = descCur;
    end

    parfor i = i_start:i_start+portion_size-1
        ii = i - i_start + 1;
        descCur = descCur_portion{ii};
        featCur = featCur_portion{ii};
        
        % now match with the specified region
        matchesModel = getMatches(descSurface, descCur, par);

        % save number of putative matches
        num_putative_all(i) = size(matchesModel, 1);

        % if number of putative matches is above a certain threshold, also save
        % the matches and the corresponding feature locations
        if num_putative_all(i) > 50
            matchesCells{i} = matchesModel; % indices of matches
            featureCells{i} = featCur; % locations of local features
            locationArray(i, :) = c; % center of large sphere
        end
    end
end
% reduce the cell arrays to contain only non-empty cells
matchesCells = matchesCells(~cellfun('isempty',matchesCells));
featureCells = featureCells(~cellfun('isempty',featureCells));
locationArray = locationArray(~cellfun('isempty',featureCells), :);
num_putative = num_putative_all(num_putative_all>50);
fprintf('Matched to all positions in %0.1f seconds...\n', toc);

%% perform RANSAC on spheres with a lot of putative matches
putative_thresh = 225;

% RANSAC options
options.minPtNum = 3; 
options.iterNum = 3e4; 
options.thDist = 0.2; 
options.thInlrRatio = 0.1; 
options.REFINE = true;
options.VERBOSE = 1;

idx = find(num_putative>putative_thresh);
num_trials = size(idx, 1);

% prepare cells to contain trial cases only
matchesCellsTrial = matchesCells(idx);
featureCellsTrial = featureCells(idx);
locationArrayTrial = locationArray(idx, :);

% little print for timing
fprintf('Performing RANSAC on all promising spheres...\n');
fprintf('This will take about %d x t minutes...\n', round(num_trials/4/60));

% statistics: numPutative, numSuccess, maxInliers, maxInlierRatio
statsPutative = zeros(num_trials, 1);
statsSuccess = zeros(num_trials, 1);
statsInliers = zeros(num_trials, 1);
statsRatio = zeros(num_trials, 1);

parfor i = 1:num_trials
    ind = idx(i);
    
    matches = matchesCellsTrial{i};
    featuresM = featureCellsTrial{i};
    pts1 = featSurface(matches(:, 1), :);
    pts2 = featuresM(matches(:, 2), :);
    
    % sanity check
    assert(num_putative(ind) == size(matches, 1));
    
    % RANSCAC
    tic
    [T, inlierPtIdx, numSuccess, maxInliers, maxInlierRatio] = ...
                ransac(pts1,pts2,options,@estimateTransform,@calcDists);
    toc
           
    % save statistics. needed to find out which sphere matches best
    
    statsPutative(i) = num_putative(ind);
    statsSuccess(i) = numSuccess;
    statsInliers(i) = maxInliers;
    statsRatio(i) = maxInlierRatio;
end

%%% EVALUTATION %%%
%% setup
centerSurface = [34.8970   21.9820   25.4546];

%% first experiment: all locations where RANSAC didn't fail
% is there a region where we have many more non-fails?
indices = find(statsSuccess > 0);
goodLocations = locationArrayTrial(indices, :);
goodLocations = [goodLocations; centerSurface];


% show all points in black and the correct center in red
color = ones(size(goodLocations, 1), 3) .* [0, 0, 0];
color(end, :) = [255, 0, 0];
pcshow(pointCloud(goodLocations, 'Color', color), 'MarkerSize', 50);


%% second experiment: more selective. which points fulfill certain
% thresholds for ALL statistics?
thSucc = 1;
thInliers = 0;
thRatio = 0; % in percent

mask = statsSuccess > thSucc & statsInliers > thInliers & statsRatio > thRatio;
mask = repmat(mask, 1, 3);
goodLocations = reshape(locationArray(mask), [], 3);
goodLocations = [goodLocations; centerSurface];

% show all points in black and the correct center in red
color = ones(size(goodLocations, 1), 3) .* [0, 0, 0];
color(end, :) = [255, 0, 0];
pcshow(pointCloud(goodLocations, 'Color', color), 'MarkerSize', 50);

%% third experiment
% look at stats of the "correct" spheres
c = centerSurface;
indices = find(statsSuccess > 0);
% get location and stats of only successful runs
statsPutativeSU = statsPutative(indices);
statsSuccessSU = statsSuccess(indices);
statsInliersSU = statsInliers(indices);
statsRatioSU = statsRatio(indices);
goodL = locationArrayTrial(indices, :);

% D is "radius" of cube
D = 4;
mask = goodL(:, 1) > c(1)-D & goodL(:, 1) < c(1)+D ...
         & goodL(:, 2) > c(2)-D & goodL(:, 2) < c(2)+D ...
         & goodL(:, 3) > c(3)-D & goodL(:, 3) < c(3)+D;     

% indices of respective locations
indices = find(mask);

mask = cat(2, mask, mask, mask);
goodLocations = reshape(goodL(mask), [], 3);    
goodLocations = [goodLocations; centerSurface];

% show all points in black and the correct center in red
color = ones(size(goodLocations, 1), 3) .* [0, 0, 0];
color(end, :) = [255, 0, 0];
pcshow(pointCloud(goodLocations, 'Color', color), 'MarkerSize', 50);

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
function mask = getDescriptorIndices(featModel, current_center, R_desc, margin)
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