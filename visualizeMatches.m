%% run GetSphericalDescriptors first
close all

% get coordinates from Model and Surface of Match
i = 67;

coordSurface = featSurface(matchesModel(inlier_idx(i), 1), :);
coordModel = featModel(matchesModel(inlier_idx(i), 2), :);

% get local points for each support region
ptsMatchSurface = getLocalPoints(ptsSurface, R, coordSurface, min_pts);
ptsMatchModel = getLocalPoints(ptsModel, R, coordModel, min_pts);

% visualize
figure('name', 'Surface');
pcshow(pointCloud(ptsMatchSurface), 'MarkerSize', 50);
figure('name', 'Model');
pcshow(pointCloud(ptsMatchModel), 'MarkerSize', 50);