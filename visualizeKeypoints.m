%% visualizeKeypoints.m
% used to create point clouds of keypoint location for presentation and for
% debugging purposes

%% load the desired model / surface
path = 'Data/PointClouds/';
pcIn = pcread(strcat(path, 'SurfaceNew_DS3.pcd'));
pts_in = pcIn.Location;

%% perform keypoint detection

% random sampling
d = 0.3;
margin = 3.5;
samples_rand = pcRandomUniformSamples(pcIn, d, margin);

% point count criterion (crit1)
min_pts = 500;
max_pts = 6000;
R = 3.5;
thVar = [3, 1.5];

mask_crit1 = false(size(samples_rand, 1), 1);
mask_crit2 = false(size(samples_rand, 1), 1);

parfor i = 1:size(samples_rand, 1)
    c = samples_rand(i, :)
    [pts_loc, dists] = getLocalPoints(pts_in, R, c, min_pts, max_pts);
    if ~isempty(dists)
        mask_crit1(i) = 1;        
        try
            [~, ~, variances] = pca(pts_loc, 'Algorithm', 'eig');
            if (variances(1) / variances(2) < thVar(1)) || ...
                    (variances(2) / variances(3) < thVar(2))
                continue
            else
                mask_crit2(i) = 1;
            end
        catch
            continue
        end        
    end
end

samples_crit1 = samples_rand(mask_crit1, :);
samples_crit2 = samples_rand(mask_crit2, :);

%% save each step as a point cloud
col_rand = [255, 0, 0].*ones(size(samples_rand, 1), 3);
col_crit1 = [0, 255, 0].*ones(size(samples_crit1, 1), 3);
col_crit2 = [0, 0, 255].*ones(size(samples_crit2, 1), 3);


pcwrite(pointCloud(samples_rand, 'Color', col_rand), 'Presentations/September 19/samples_rand.pcd');
pcwrite(pointCloud(samples_crit1, 'Color', col_crit1), 'Presentations/September 19/samples_crit1.pcd');
pcwrite(pointCloud(samples_crit2, 'Color', col_crit2), 'Presentations/September 19/samples_crit2.pcd');

%% save each stage in combination with Surface Point Cloud
pts_comb_rand = [pts_in; samples_rand];
pc_comb_rand = pointCloud(pts_comb_rand, 'Color', [pcIn.Color; col_rand]);
pcwrite(pc_comb_rand, 'Presentations/September 19/samples_rand_comb.pcd');

pts_comb_crit1 = [pts_in; samples_crit1];
pc_comb_crit1 = pointCloud(pts_comb_crit1, 'Color', [pcIn.Color; col_crit1]);
pcwrite(pc_comb_crit1, 'Presentations/September 19/samples_crit1_comb.pcd');

pts_comb_crit2 = [pts_in; samples_crit2];
pc_comb_crit2 = pointCloud(pts_comb_crit2, 'Color', [pcIn.Color; col_crit2]);
pcwrite(pc_comb_crit2, 'Presentations/September 19/samples_crit2_comb.pcd');

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