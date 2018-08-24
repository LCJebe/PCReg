%% run getSphericalDescriptors first to obtain putative matches

addpath('../');



%% Experiment 2: After matching: Putative matches & inlier matches
% This experiemnt is meant to get statistics about the output of the 
% matching algorithm By visualizing properties such as # points, thVar and
% avgDist from center, hopefully a more informed decision about keypoint
% detection / sampling parameters can be made. 

INLIER = false;

point_countsS = [];
avg_distsS = [];
variance_ratio1S = [];
variance_ratio2S = [];
point_countsM = [];
avg_distsM = [];
variance_ratio1M = [];
variance_ratio2M = [];

if INLIER
    matches = matchesModel(inlier_idx, :);
else
    matches = matchesModel;
end

for i = 1:size(matches, 1)
    % matches
    mS = matchesModel(i, 1);
    mM = matchesModel(i, 2);
    
    % features 
    fS = featSurface(mS, :);
    fM = featModel(mM, :);
    
    % given the feature location, get the local points that were used for this feature
    % points and distances
    [pS, dS] = getLocalPoints(pcSurface.Location, descOpt.R, fS, 0, inf);
    [pM, dM] = getLocalPoints(pcModel.Location, descOpt.R, fM, 0, inf);
    
    %%%%%%%%%%%%%%%
    % collect statistics: Surface
    %%%%%%%%%%%%%%%
    
    % 1) get variance ratios
    [~, ~, var] = pca(pS, 'Algorithm', 'eig');
    var1 = var(1) / var(2);
    var2 = var(2) / var(3);
    variance_ratio1S = [variance_ratio1S, var1];
    variance_ratio2S = [variance_ratio2S, var2];
    
    % 2) number of points
    point_countsS = [point_countsS, size(pS, 1)];
    
    % 3) average distance from center
    avg_distsS = [avg_distsS, mean(dS)];
    
    %%%%%%%%%%%%%%%
    % collect statistics: Model
    %%%%%%%%%%%%%%%
    
    % 1) get variance ratios
    [~, ~, var] = pca(pM, 'Algorithm', 'eig');
    var1 = var(1) / var(2);
    var2 = var(2) / var(3);
    variance_ratio1M = [variance_ratio1M, var1];
    variance_ratio2M = [variance_ratio2M, var2];
    
    % 2) number of points
    point_countsM = [point_countsM, size(pM, 1)];
    
    % 3) average distance from center
    avg_distsM = [avg_distsM, mean(dM)];
end
   
   
%% Plot histograms
close all

% window options
% figure settings
screensize = get( 0, 'Screensize' );
figpos = [screensize(3)/6, 125, 2*screensize(3)/3, screensize(4)-250];


binlim_counts = [0, 10000];
binlim_dists = [1.5, 3.5];
binlim_var1 = [1, 10];
binlim_var2 = [1, 10];


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