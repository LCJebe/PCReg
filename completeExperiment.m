%% slideMatchingWindow.m
% experiment, where we slowly slide the crop from the ideal position
% outwards with a ramdom angle, and record how well RANSAC performs with
% different metrics (inlier percentage, num successes, num inliers)

rng('shuffle');

%% read Model and Surface pointcloud
path = 'Data/PointClouds/';
pcModel = pcread(strcat(path, 'ModelSmoothColorUp3.pcd'));
pcSurface = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));

% shift pointclouds to center
pcSurface = centerPointCloud(pcSurface);
pcModel = centerPointCloud(pcModel);

%% get descriptors for surface and for the complete model
dS = 0.3; % denser (0.2 or 0.3)
dM = 0.4; % sparser (0.4 or more)
margin = 3.5;

sample_ptsSurface = pcRandomUniformSamples(pcSurface, dS, margin);
sample_ptsModel = pcRandomUniformSamples(pcModel, dM, -margin);

% descriptor options
descOptM.ALIGN_POINTS = true;
descOptM.CENTER = false;
descOptM.min_pts = 500;
descOptM.max_pts = 6000;
descOptM.R = 3.5;
descOptM.thVar = [3, 1.5]; 
descOptM.k = 0.85;
descOptM.VERBOSE = 1;
descOptM.max_region_size = 10;

descOptS = descOptM;
descOptS.max_region_size = 100;

[featModel, descModel] = ...
        speedyDescriptors(pcModel.Location, sample_ptsModel, descOptM);   
% save the features and descriptors to workspace
save('featModel0.4.mat', 'featModel');
save('descModel0.4.mat', 'descModel');

[featSurface, descSurface] = ...
        getSpacialHistogramDescriptors(pcSurface.Location, sample_ptsSurface, descOptS);   

    
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

tic
parpool()
fprintf('Loaded parpool in %0.1f seconds...\n', toc);

% first: check for valid spheres by number of descriptors in sphere
tic
min_pts = 50;
max_pts = inf;
valid_mask = false(size(sample_ptsSpheres, 1), 1);
parfor i = 1:size(sample_ptsSpheres, 1)
    c = sample_ptsSpheres(i, :);
    [~, dists] = getLocalPoints(featModel, R_desc, c, min_pts, max_pts);
    if ~isempty(dists)
        valid_mask(i) = 1;
    end
end
fprintf('Obtained test positions in %0.1f seconds...\n', toc);

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

num_putative = zeros(size(sample_ptsSpheres, 1), 1);
matchesCells = cell(size(sample_ptsSpheres, 1), 1);
featureCells = cell(size(sample_ptsSpheres, 1), 1);
locationCells = cell(size(sample_ptsSpheres, 1), 1);

tic
num_spheres = size(sample_ptsSpheres, 1);
fprintf('This will take about %d minutes...\n', round(num_spheres/4/60));
parfor i = 1:num_spheres
    f = waitbar_handle;
    c = sample_ptsSpheres(i, :);
    
    % mask of descriptors that lie within sphere
    mask = getDescriptorIndices(featModel, c, R_desc, 0);
        
    % based on the mask, return the relevant sample points and descriptors
    maskF = cat(2, mask, mask, mask);
    featCur = reshape(featModel(maskF), [], 3);  
    maskD = repmat(mask, 1, size(descModel, 2));
    descCur = reshape(descModel(maskD), [], size(descModel, 2));
    
    % now match with the specified region
    matchesModel = getMatches(descSurface, descCur, par);
    
    % save number of putative matches
    num_putative(i) = size(matchesModel, 1);
    
    % if number of putative matches is above a certain threshold, also save
    % the matches and the corresponding feature locations
    if num_putative(i) > 50
        matchesCells{i} = matchesModel;
        featureCells{i} = featCur;
        locationCells{i} = c;
    end
    
end

% reduce the cell arrays to contain only non-empty cells
matchesCells = matchesCells(~cellfun('isempty',matchesCells));
featureCells = featureCells(~cellfun('isempty',featureCells));
locationCells = locationCells(~cellfun('isempty',locationCells));
fprintf('Matches to all positions in %0.1f seconds...\n', toc);

%% perform RANSAC on spheres with a lot of putative matches


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

%% helper function: descriptors within certain distance of given center point
function mask = getDescriptorIndices(featModel, current_center, R_desc, margin)
    feat_rel = featModel - current_center;
    dists = vecnorm(feat_rel, 2, 2);
    mask = dists < (R_desc + margin); % note that margin < 0
end