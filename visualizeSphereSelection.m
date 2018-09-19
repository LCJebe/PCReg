%% load complete experiment workspace to start
clear all
load('completeExperimentFast_workspace_d5m0.3s0.3_raw.mat');


%%
% create gray pointcloud only
col1 = 0.4 * ones(size(sample_ptsSpheres, 1), 3);
pc1 = pointCloud(single(sample_ptsSpheres), 'Color', col1);

% positions with most putative matches: locationArrayTrial
putative_thresh = 250;
idx = find(num_putativeTrial>putative_thresh);
locationArrayTrialNew = locationArrayTrial(idx, :);

statsPutativeNew = statsPutative(idx);
statsSuccessNew = statsSuccess(idx);
statsInliersNew = statsInliers(idx);
statsRatioNew = statsRatio(idx);

% create gray and blue pointcloud
col2 = [0, 0, 1] .* ones(size(locationArrayTrialNew, 1), 3);
col12 = [col1; col2];
loc12 = [sample_ptsSpheres; locationArrayTrialNew];
pc12 = pointCloud(single(loc12), 'Color', col12);

% positions where RANSAC Succeded
thSucc = 0;
thInliers = 28;
thRatio = 10; % in percent
thPutative = 170;

mask = statsSuccess >= thSucc & statsInliers >= thInliers & statsRatio >= thRatio & statsPutative >= thPutative;
indices = find(mask);
goodLocations = locationArrayTrial(indices, :);

% create gray, blue, and red pointcloud
col3 = [1, 0, 0] .* ones(size(goodLocations, 1), 3);
loc123 = [loc12; goodLocations];
col123 = [col12; col3];
pc123 = pointCloud(single(loc123), 'Color', col123);

% create red pointcloud only
pc3 = pointCloud(single(goodLocations), 'Color', col3);
size(col3)

%pcshow(pc3, 'MarkerSize', 50);

%% save all pointclouds to presentation folder
path = 'Presentations/September 19/';

%pcwrite(pc1, strcat(path, 'fastpc1.pcd'));
%pcwrite(pc12, strcat(path, 'fastpc12.pcd'));
%pcwrite(pc123, strcat(path, 'fastpc123.pcd'));
%pcwrite(pc3, strcat(path, 'fastpc3.pcd'));

