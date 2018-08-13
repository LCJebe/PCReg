%% This script performs a series of tests that investigate the 
% RANSAC / MATCHING / LOCAL ALIGNMENT Problem
clear all
close all

%% Experiment 1) check rotation invariance of descriptors directly

% get the supporting points for one single descriptor
pc = pcread('teapot.ply');

pts = pc.Location;
%c = [3.3, 0, 2.4]; % support region center
c = [2, -0.5, 1]; % support region center
R = 1; % support radius

pts_rel = pts - c; % relative to center
dists = vecnorm(pts_rel, 2, 2);
mask = dists < R;
dists = dists(find(mask));
mask = cat(2, mask, mask, mask);
pts_sphere = reshape(pts_rel(mask), [], 3);    

% now create N random rotations of the points
D = 31;
N = 30;
points = cell(N, 1);

for i = 1:N
    rotation = rand(1, 3)*2*pi;
    Rotm = eul2rotm(rotation);
    pts_rot = pts_sphere*Rotm;
    points{i} = pts_rot;
end

min_pts = 30;
sample_pts = [0, 0, 0]; % points are already local
thVar = [1, 1];
ALIGN = true;
CENTER = false;

% get descriptors
descriptors = zeros(N, D);
for i = 1:N
    pts_rot = points{i};
    [feat, desc] = getMomentDescriptors(pts_rot, sample_pts, min_pts, R, thVar, ALIGN, CENTER, 'all');
    descriptors(i, :) = desc;
end


%% Experiment 2: Match the Surface with a transformed version of itself
path = 'Data/PointClouds/';
pcSurface = pcread(strcat(path, 'Surface_DS2.pcd'));

% specify rotation and translation
r = [0.4, 1, -1.2];
t = [-3, 8, 25];
Rotm = eul2rotm(r);

T = eye(4);
T(1:3, 1:3) = Rotm;
T(4, 1:3) = t;

pcSurface_rot = pctransform(pcSurface, affine3d(T));

% center both pointclouds
%pcSurface = centerPointCloud(pcSurface);
%pcSurface_rot = centerPointCloud(pcSurface_rot);

% get sample points
d = 0.4;
sample_ptsSurface = pcRandomUniformSamples(pcSurface, d);
sample_ptsSurface_rot = pcRandomUniformSamples(pcSurface_rot, d);

% define parameters and get descriptors
ptsSurface = pcSurface.Location;
ptsSurface_rot = pcSurface_rot.Location;
min_pts = 50;
R = 1.5;
thVar = [1.5, 1.5];
ALIGN = true;
CENTER = false;
k = 'all';

[featSurface, descSurface] = getMomentDescriptors(ptsSurface, sample_ptsSurface, min_pts, R, thVar, ALIGN, CENTER, k);
[featSurface_rot, descSurface_rot] = getMomentDescriptors(ptsSurface_rot, sample_ptsSurface_rot, min_pts, R, thVar, ALIGN, CENTER, k);

% apply weights
weights.M1 = 0; % 0
weights.M2 = 0.3; % 0.3
weights.M3 = 0.8; % 0.8
weights.M4 = 0.5; % 0.5

descSurface = applyMomentWeight(descSurface, weights);
descSurface_rot = applyMomentWeight(descSurface_rot, weights);

% define matchin parameters and match descriptors

par.Method = 'Approximate'; % 'Exhaustive' (default) or 'Approximate'
par.MatchThreshold =  0.55; % 1.0 (default) Percent Value (0 - 100) for distance-reject
par.MaxRatio = 0.6; % 0.6 (default) nearest neighbor ambiguity rejection
par.Metric =  'SSD'; % SSD (default) for L2, SAD for L1
par.Unique = true; % true: 1-to-1 mapping only, else set false (default)

matches = matchFeatures(descSurface, descSurface_rot, ...
        'Method', par.Method, ...
        'MatchThreshold', par.MatchThreshold, ... 
        'MaxRatio', par.MaxRatio, ... 
        'Metric', par.Metric, ...
        'Unique', par.Unique); 
    
% get location of matches points
loc1M = featSurface(matches(:, 1), :);
loc1S = featSurface_rot(matches(:, 2), :);

% now call getInliersRANSAC to check if the alignment succeeded

%% now use T from RANSAC and check, if the pointclouds are actually aligned
pts = pcSurface.Location;
color = pcSurface.Color;
pts_trans = [pts, ones(size(pts, 1), 1)] * f1;
pcSurface_trans = pointCloud(pts_trans(:, 1:3), 'Color', color);
figure()
pcshow(pcSurface_trans);
figure()
pcshow(pcSurface_rot);


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