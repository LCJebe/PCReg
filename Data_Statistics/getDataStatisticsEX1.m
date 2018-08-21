addpath('../');

%% load model and surface
path = '../Data/PointClouds/';

pcSurface = pcread(strcat(path, 'Surface_DS3_alignedM.pcd'));
pcModel = pcread(strcat(path, 'GoodCropSmoothUp3_large.pcd'));

%% Experiment 1: Before matching
% This experiemnt is meant to get statistics about the input to the
% matching algorithm (or rather, statistics even before descriptor
% calculation). By visualizing properties such as # points, thVar and
% avgDist from center, hopefully a more informed decision about keypoint
% detection / sampling parameters can be made. 

% keypoint detector options
d = 1;
marginSurface = 3.5;
marginModel = -3.5;
sample_ptsSurface = pcRandomUniformSamples(pcSurface, d, marginSurface);
sample_ptsModel = pcRandomUniformSamples(pcModel, d, marginModel);

R = 3.5;
min_pts = 200;
max_pts = inf;
thVar = [1, 1];

% Surface
point_countsS = [];
avg_distsS = [];
variance_ratio1S = [];
variance_ratio2S = [];
tic
for i = 1:size(sample_ptsSurface ,1)
    c = sample_ptsSurface(i, :);
    [local_pts, dists] = getLocalPoints(pcSurface.Location, R, c, min_pts, max_pts);
    if ~isempty(local_pts)        
        % get variance ratios
        [~, ~, var] = pca(local_pts, 'Algorithm', 'eig');
        var1 = var(1) / var(2);
        var2 = var(2) / var(3);
        if var1 >= thVar(1) && var2 >= thVar(2)
            % save number of points
            point_countsS = [point_countsS, size(local_pts, 1)];
            % save average distance from center
            avg_distsS = [avg_distsS, mean(dists)];
            % save variance ratio
            variance_ratio1S = [variance_ratio1S, var1];
            variance_ratio2S = [variance_ratio2S, var2];
        end
    end
end

% little information print
num_valid = length(point_countsS);
fprintf('Collected statistics for %d valid points (%0.2f %%) in %0.2f seconds for Surface\n', ...
       num_valid, 100*num_valid/size(sample_ptsSurface ,1), toc);


% Model
point_countsM = [];
avg_distsM = [];
variance_ratio1M = [];
variance_ratio2M = [];
tic
for i = 1:size(sample_ptsModel ,1)
    c = sample_ptsModel(i, :);
    [local_pts, dists] = getLocalPoints(pcModel.Location, R, c, min_pts, max_pts);
    if ~isempty(local_pts)
        % get variance ratios
        [~, ~, var] = pca(local_pts, 'Algorithm', 'eig');
        var1 = var(1) / var(2);
        var2 = var(2) / var(3);
        if var1 >= thVar(1) && var2 >= thVar(2)
            % save number of points
            point_countsM = [point_countsM, size(local_pts, 1)];
            % save average distance from center
            avg_distsM = [avg_distsM, mean(dists)];
            % save variance ratio
            variance_ratio1M = [variance_ratio1M, var1];
            variance_ratio2M = [variance_ratio2M, var2];
        end
    end
end

% little information print
num_valid = length(point_countsM);
fprintf('Collected statistics for %d valid points (%0.2f %%) in %0.2f seconds for Model\n', ...
       num_valid, 100*num_valid/size(sample_ptsModel ,1), toc);
   
   
%% Plot histograms
close all

% window options
% figure settings
screensize = get( 0, 'Screensize' );
figpos = [screensize(3)/6, 125, 2*screensize(3)/3, screensize(4)-250];

BIVARIATE = false;


% define bin limits that are valid for both Surface and Model --> easier
% comparison

binlim_counts = [0, 1.5e4];
binlim_dists = [1.5, 3.5];
binlim_var1 = [0, 5];
binlim_var2 = [0, 5];


% Surface
fig_h = figure();
set(fig_h,'Position',figpos)
subplot(2, 2, 1);
histogram(point_countsS, 'BinLimits', binlim_counts);
title('Surface: Number of points in local neighborhood');
ylabel('Occurrences');
xlabel('# points');

subplot(2, 2, 2);
histogram(avg_distsS, 'BinLimits', binlim_dists);
title('Surface: Average distance of points to center within one local neighborhood');
ylabel('Occurrences');
xlabel('avg. distance');

if BIVARIATE
    edges = {0:0.5:5, 0:0.5:5};
    subplot(2, 1, 2);
    hist3([variance_ratio1S', variance_ratio2S'], 'Edges', edges);
    title('Surface: Variance Ratio 1');
    ylabel('Occurrences');
    xlabel('Variance Ratio 1');
else
    subplot(2, 2, 3);
    histogram(variance_ratio1S, 'BinLimits', binlim_var1);
    title('Surface: Variance Ratio 1');
    ylabel('Occurrences');
    xlabel('Variance Ratio 1');

    subplot(2, 2, 4);
    histogram(variance_ratio2S, 'BinLimits', binlim_var2);
    title('Surface: Variance Ratio 2');
    ylabel('Occurrences');
    xlabel('Variance Ratio 2');
end

% Model
fig_h = figure();
set(fig_h,'Position',figpos)
subplot(2, 2, 1);
histogram(point_countsM, 'BinLimits', binlim_counts);
title('Model: Number of points in local neighborhood');
ylabel('Occurrences');
xlabel('# points');

subplot(2, 2, 2);
histogram(avg_distsM, 'BinLimits', binlim_dists);
title('Model: Average distance of points to center within one local neighborhood');
ylabel('Occurrences');
xlabel('avg. distance');

subplot(2, 2, 3);
histogram(variance_ratio1M, 'BinLimits', binlim_var1);
title('Model: Variance Ratio 1');
ylabel('Occurrences');
xlabel('Variance Ratio 1');

subplot(2, 2, 4);
histogram(variance_ratio2M, 'BinLimits', binlim_var2);
title('Model: Variance Ratio 2');
ylabel('Occurrences');
xlabel('Variance Ratio 2');


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