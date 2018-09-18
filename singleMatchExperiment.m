%% load model and surface descriptors with locations

load('Data/Descriptors/featSurface0.3_raw.mat');
load('Data/Descriptors/descSurface0.3_raw.mat'); 
load('Data/Descriptors/featModel0.3.mat');
load('Data/Descriptors/descModel0.3.mat'); 

%% match only once (uniquely), and check distribution of matches
% matching options
par.UNNORMALIZE = true;
par.norm_factor = 2; % 2

par.CHANGE_METRIC = true;
par.metric_factor = 0.6; % 0.45 / 0.6

par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold = 1; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.60; % 0.6 (default) nearest neighbor ambiguity rejection
par.Metric =  'SAD'; % SSD (default) for L2, SAD for L1
par.Unique = true; % true: 1-to-1 mapping only, else set false (default)

par.VERBOSE = 1;

tic
matchesModel = getMatches(descSurface, descModel, par);
toc

%% show matched positions of model in 3D and center Surface in red
centerSurface = [34.8970   21.9820   25.4546];
 
matchesPoints = [featModel(matchesModel(:, 2), :); centerSurface];

col = [80, 80, 80] .* ones(size(matchesPoints, 1), 3) / 255;
col(end, :) = [1, 0, 0];

pcshow(pointCloud(matchesPoints, 'Color', col), 'Markersize', 50);