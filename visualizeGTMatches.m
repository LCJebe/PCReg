%% visualize GT matches
% this code takes ground truth correspondeces from the good crop and
% aligned Model, and tries to align them using PCA, Sign disambiguation,
% and optional %KNN or local points only

close all

%% load model and aligned surface
path = 'Data/PointClouds/';

pcSurface = pcread(strcat(path, 'SurfaceNew_DS3_alignedM.pcd'));
pcModel = pcread(strcat(path, 'GoodCropSpherical_smaller.pcd'));

%% options
MAX_MATCHES = 0;
METHOD = 'KNN'; % 'C', 'KNN', KNN_C' % just for the display of points
TEST_KNN = false; % compare different KNN methods (for Sept 19 Presentation)
TEST_C = false; % test effect of center vs. centroid in alignment

%% sample points (only intersection of surface and model sampling area)
d = 0.3;
margin = 3.5;
sample_pts = pcIntersectionSamples(pcSurface, pcModel, d, margin);

%% select points, where surface and model are not empty
R = 3.5;
min_pts = 500;
max_pts = 6000;
thVar = [3, 1.5];

%% determine the valid points and sample only those
% only when both sopport regions have Nmin < Npoints < Nmax points
mask_crit1 = false(size(sample_pts, 1), 1);
mask_crit2 = false(size(sample_pts, 1), 1);
ptsSurface = pcSurface.Location;
ptsModel = pcModel.Location;

tic
local_ptsSCell = cell(size(sample_pts, 1), 1);
local_ptsMCell = cell(size(sample_pts, 1), 1);
ratiosSArray = zeros(size(sample_pts, 1), 2);
ratiosMArray = zeros(size(sample_pts, 1), 2);
locationArray = zeros(size(sample_pts, 1), 3);

parfor i = 1:size(sample_pts, 1)
    c = sample_pts(i, :);
    [local_ptsS, distsS] = getLocalPoints(ptsSurface, R, c, min_pts, max_pts);
    [local_ptsM, distsM] = getLocalPoints(ptsModel, R, c, min_pts, max_pts);
    if ~isempty(distsS) && ~isempty(distsM)
        mask_crit1(i) = 1;
        try            
            [~, ~, varsS] = pca(local_ptsS, 'Algorithm', 'eig');        
            [~, ~, varsM] = pca(local_ptsM, 'Algorithm', 'eig');
            ratiosS = [varsS(1) / varsS(2), varsS(2) / varsS(3)];
            ratiosM = [varsM(1) / varsM(2), varsM(2) / varsM(3)];

            if all((ratiosS - thVar) > 0) && all((ratiosM - thVar) > 0)
                mask_crit2(i) = 1;
                local_ptsSCell{i} = local_ptsS;
                local_ptsMCell{i} = local_ptsM;
                ratiosSArray(i, :) = ratiosS;
                ratiosMArray(i, :) = ratiosM;
                locationArray(i, :) = c;
            else
                continue;
            end
        catch
            continue
        end      
    end
end

% apply mask
valid_pts = sample_pts(mask_crit2, :);
local_ptsSCell = local_ptsSCell(mask_crit2);
local_ptsMCell = local_ptsMCell(mask_crit2);
ratiosSArray = ratiosSArray(mask_crit2, :);
ratiosMArray = ratiosMArray(mask_crit2, :);
locationArray = locationArray(mask_crit2, :);

valid_numpts = [cellfun(@(a) size(a, 1), local_ptsSCell), cellfun(@(a) size(a, 1), local_ptsMCell)];
valid_vars = [ratiosSArray, ratiosMArray];
desc_dists = [];
desc_dists2 = [];
desc_dists3 = [];
desc_dists4 = [];
alignment_meas = [];
alignment_meas2 = [];
alignment_meas3 = [];
alignment_meas4 = [];
centroid_dists2 = [];

%%
for i = 1:size(valid_pts ,1)
    
    % unpack local points
    local_ptsS = local_ptsSCell{i};
    local_ptsM = local_ptsMCell{i};
    
    if ~TEST_KNN && ~ TEST_C
        % align using both methods and save descriptor distance
        [local_ptsS_aligned, tfS] = AlignPoints(local_ptsS);
        [local_ptsM_aligned, tfM] = AlignPoints(local_ptsM);

        [local_ptsS_aligned2, tfS2, cS2] = AlignPoints_c(local_ptsS);
        [local_ptsM_aligned2, tfM2, cM2] = AlignPoints_c(local_ptsM);

        % DEBUG!!! use captialized KNN not lower case knn!!
        [local_ptsS_aligned3, tfS3, cS3] = AlignPoints_KNN(local_ptsS);
        [local_ptsM_aligned3, tfM3, cM3] = AlignPoints_KNN(local_ptsM);    

        [local_ptsS_aligned4, tfS4, cS4] = AlignPoints_KNN_c(local_ptsS);
        [local_ptsM_aligned4, tfM4, cM4] = AlignPoints_KNN_c(local_ptsM);
    elseif TEST_KNN
        % align using both methods and save descriptor distance
        [local_ptsS_aligned, tfS] = AlignPoints_KNN(local_ptsS);
        [local_ptsM_aligned, tfM] = AlignPoints_KNN(local_ptsM);

        k = 500;
        [local_ptsS_aligned2, tfS2, cS2] = AlignPoints_knn(local_ptsS, k);
        [local_ptsM_aligned2, tfM2, cM2] = AlignPoints_knn(local_ptsM, k);
 
        k = 1500;
        [local_ptsS_aligned3, tfS3, cS3] = AlignPoints_knn(local_ptsS, k);
        [local_ptsM_aligned3, tfM3, cM3] = AlignPoints_knn(local_ptsM, k);

        [local_ptsS_aligned4, tfS4, cS4] = AlignPoints_weighted(local_ptsS);
        [local_ptsM_aligned4, tfM4, cM4] = AlignPoints_weighted(local_ptsM);
    elseif TEST_C
        % align using both methods and save descriptor distance
        [local_ptsS_aligned, tfS] = AlignPoints_KNN(local_ptsS, false, false);
        [local_ptsM_aligned, tfM] = AlignPoints_KNN(local_ptsM, false, false);
        
        [local_ptsS_aligned2, tfS2] = AlignPoints_KNN(local_ptsS, true, false);
        [local_ptsM_aligned2, tfM2] = AlignPoints_KNN(local_ptsM, true, false);
        
        [local_ptsS_aligned3, tfS3] = AlignPoints_KNN(local_ptsS, false, true);
        [local_ptsM_aligned3, tfM3] = AlignPoints_KNN(local_ptsM, false, true);
        
        [local_ptsS_aligned4, tfS4] = AlignPoints_KNN(local_ptsS, true, true);
        [local_ptsM_aligned4, tfM4] = AlignPoints_KNN(local_ptsM, true, true);
    end

    if ~isempty(local_ptsS_aligned2) && ~isempty(local_ptsM_aligned2)

        % get descriptors for both alignments
        dS = getDesc(local_ptsS_aligned, R);
        dM = getDesc(local_ptsM_aligned, R);

        dS2 = getDesc(local_ptsS_aligned2, R);
        dM2 = getDesc(local_ptsM_aligned2, R);

        dS3 = getDesc(local_ptsS_aligned3, R);
        dM3 = getDesc(local_ptsM_aligned3, R);
        
        dS4 = getDesc(local_ptsS_aligned4, R);
        dM4 = getDesc(local_ptsM_aligned4, R);

        % get descriptor distance for both alignments
        dist = getDescDist(dS, dM);
        dist2 = getDescDist(dS2, dM2);
        dist3 = getDescDist(dS3, dM3);
        dist4 = getDescDist(dS4, dM4);

        % save in complete variables
        desc_dists = [desc_dists, dist];
        desc_dists2 = [desc_dists2, dist2];
        desc_dists3 = [desc_dists3, dist3];
        desc_dists4 = [desc_dists4, dist4];

        % get measurement for how good the alignment is
        alignment_meas = [alignment_meas, checkAlignment(tfS, tfM)];
        alignment_meas2 = [alignment_meas2, checkAlignment(tfS2, tfM2)];
        alignment_meas3 = [alignment_meas3, checkAlignment(tfS3, tfM3)];
        alignment_meas4 = [alignment_meas4, checkAlignment(tfS4, tfM4)];

        % save centroid distance for alignment 2
        centroid_dists2 = [centroid_dists2, norm(cS2-cM2)];
    else
        % alignment_v2 failed, but we still need to fill the
        % vectors with dummies so that indexing doesn't break
        % "nan" can be ignored with nanmedian and nanmean later!!
        desc_dists = [desc_dists, nan];
        desc_dists2 = [desc_dists2, nan];
        desc_dists3 = [desc_dists3, nan];
        desc_dists4 = [desc_dists4, nan];
        alignment_meas = [alignment_meas, nan];
        alignment_meas2 = [alignment_meas2, nan];
        alignment_meas3 = [alignment_meas3, nan];
        alignment_meas4 = [alignment_meas4, nan];
        
        centroid_dists2 = [centroid_dists2, nan];
    end
end
toc
clear local_ptsS local_ptsM

% print out some intermediate statistics: valid samples points based on
% constraints
num_samples = size(sample_pts ,1);
num_valid = size(valid_pts, 1);
perc_valid = double(num_valid) / num_samples * 100;
fprintf('Sampled %d points, %d are valid (%0.2f %%)\n', ...
   num_samples, num_valid, perc_valid);

% and some more cool stuff: descriptor distances
median_dist = nanmedian(desc_dists);
median_dist2 = nanmedian(desc_dists2);
median_dist3 = nanmedian(desc_dists3);
median_dist4 = nanmedian(desc_dists4);
mean_dist = nanmean(desc_dists);
mean_dist2 = nanmean(desc_dists2);
mean_dist3 = nanmean(desc_dists3);
mean_dist4 = nanmean(desc_dists4);
fprintf('Median (Mean) Distance with old alignment: %0.2f (%0.2f)\n', median_dist, mean_dist);
fprintf('Median (Mean) Distance with local alignment: %0.2f (%0.2f)\n', median_dist2, mean_dist2);
fprintf('Median (Mean) Distance with knn alignment: %0.2f (%0.2f)\n', median_dist3, mean_dist3);
fprintf('Median (Mean) Distance with knn and local alignment: %0.2f (%0.2f)\n', median_dist4, mean_dist4);

%% and more cool statistics: num successes 
successThresh = 0.5;
num_success = sum(alignment_meas < successThresh);
num_success2 = sum(alignment_meas2 < successThresh);
num_success3 = sum(alignment_meas3 < successThresh);
num_success4 = sum(alignment_meas4 < successThresh);
perc_success = 100*num_success/num_valid;
perc_success2 = 100*num_success2/num_valid;
perc_success3 = 100*num_success3/num_valid;
perc_success4 = 100*num_success4/num_valid;
fprintf('Old: %d (%0.1f %%) successes, Local: %d (%0.1f %%) successes, KNN: %d (%0.1f %%) successes\n', ...
            num_success, perc_success, num_success2, perc_success2, num_success3, perc_success3);
        
PLOT = true;
if PLOT
    thresh = 0:0.001:0.5;
    perc_succ = zeros(length(thresh), 1);
    perc_succ2 = zeros(length(thresh), 1);
    perc_succ3 = zeros(length(thresh), 1);
    perc_succ4 = zeros(length(thresh), 1);
    for i = 1:length(thresh)
        th = thresh(i);
        num_success = sum(alignment_meas < th);
        num_success2 = sum(alignment_meas2 < th);
        num_success3 = sum(alignment_meas3 < th);
        num_success4 = sum(alignment_meas4 < th);
        perc_succ(i) = 100*num_success/num_valid;
        perc_succ2(i) = 100*num_success2/num_valid;
        perc_succ3(i) = 100*num_success3/num_valid;
        perc_succ4(i) = 100*num_success4/num_valid;
    end
    plot(thresh, perc_succ);
    hold all
    plot(thresh, perc_succ3);
    plot(thresh, perc_succ2);
    plot(thresh, perc_succ4);
    if ~TEST_KNN && ~TEST_C
        legend('Global PCA', 'KNN PCA, 0.85%', 'Local PCA, r = 2.5', 'KNN & Local PCA', 'Location', 'NorthWest');
    elseif TEST_KNN
        legend('KNN PCA, 0.85%', 'k = 500', 'k = 1500', 'weighted points', 'Location', 'NorthWest');
    elseif TEST_C
        legend('KNN PCA, 0.85%', 'PCA center, SGN centroid', 'SGN center, PCA centroid', 'PCA & SGN from center', 'Location', 'NorthWest');
    end
    xlabel('threshold');
    ylabel('sucess rate');
    grid;
end
        
%% more stuff: centroid recognition distance
fprintf('Median (Mean) Centroid Distance: %0.2f, (%0.2f)\n', ...
            nanmedian(centroid_dists2), nanmean(centroid_dists2));


%% now the fun part: visualize ground truth matches

% for better display: generate the 8 points that span the cube defined by R
add_pts = R*[-1,-1,-1;...
             -1,-1, 1;...
             -1, 1,-1;...
             -1, 1, 1;...
              1,-1,-1;...
              1,-1, 1;...
              1, 1,-1;...
              1, 1, 1];

% figure settings
screensize = get( 0, 'Screensize' );
figpos = [screensize(3)/8, 75, 3*screensize(3)/4, screensize(4)-150];

num_display = min([MAX_MATCHES, num_valid]);

% get a random selection of the matches
matches_vis = randperm(size(valid_pts, 1),num_display);

for i = 1:num_display
    idx = matches_vis(i);

    c = valid_pts(idx, :);
    
    % get points
    [ptsMatchSurface, ~] = getLocalPoints(pcSurface.Location, R, c, min_pts, max_pts);
    [ptsMatchModel, ~] = getLocalPoints(pcModel.Location, R, c, min_pts, max_pts);
    
    % align to local reference frame using both methods       
    [ptsMatchSurface_aligned, ~] = AlignPoints(ptsMatchSurface);
    [ptsMatchModel_aligned, ~] = AlignPoints(ptsMatchModel);

    % perform the method specified for visual comparison
    if strcmp(METHOD, 'C')
        [ptsMatchSurface_aligned_disp, ~, cS] = AlignPoints_c(ptsMatchSurface);
        [ptsMatchModel_aligned_disp, ~, cM] = AlignPoints_c(ptsMatchModel);
    elseif strcmp(METHOD, 'KNN')
        [ptsMatchSurface_aligned_disp, ~, cS] = AlignPoints_KNN(ptsMatchSurface);
        [ptsMatchModel_aligned_disp, ~, cM] = AlignPoints_KNN(ptsMatchModel);
    elseif strcmp(METHOD, 'KNN_C')
        [ptsMatchSurface_aligned_disp, ~, cS] = AlignPoints_KNN_c(ptsMatchSurface);
        [ptsMatchModel_aligned_disp, ~, cM] = AlignPoints_KNN_c(ptsMatchModel);
    end

    
    % add frame points for display
    ptsMatchSurface = [ptsMatchSurface; add_pts];
    ptsMatchModel = [ptsMatchModel; add_pts];
    ptsMatchSurface_aligned = [ptsMatchSurface_aligned; add_pts];
    ptsMatchModel_aligned = [ptsMatchModel_aligned; add_pts];
    ptsMatchSurface_aligned_disp = [ptsMatchSurface_aligned_disp; add_pts];
    ptsMatchModel_aligned_disp = [ptsMatchModel_aligned_disp; add_pts];
    
    % make plots
    row = mod(i-1, 4)+1;
    if row == 1
       fig_h = figure();
       set(fig_h,'Position',figpos)
    end
    subplot(4, 6, 6*row-5);
    pcshow(ptsMatchSurface, 'MarkerSize', 50);
    xlabel(sprintf('%d pts', valid_numpts(idx, 1)));
    ylabel(sprintf('cd: %0.2f', norm(cS-cM))); % centroid distance
    if row == 1; title('Surface'); end
    subplot(4, 6, 6*row-4);
    pcshow(ptsMatchModel, 'MarkerSize', 50);
    xlabel(sprintf('%d pts', valid_numpts(idx, 2)));
    if row == 1; title('Model'); end
    subplot(4, 6, 6*row-3);
    pcshow(ptsMatchSurface_aligned, 'MarkerSize', 50);
    xlabel(sprintf('align: %0.2f', alignment_meas(idx)));
    if row == 1; title('Surface aligned'); end
    subplot(4, 6, 6*row-2);
    pcshow(ptsMatchModel_aligned, 'MarkerSize', 50);
    if row == 1; title('Model aligned'); end
    subplot(4, 6, 6*row-1);
    pcshow(ptsMatchSurface_aligned_disp, 'MarkerSize', 50);
    xlabel(sprintf('align: %0.2f', alignment_meas2(idx)));
    if row == 1; title(strcat('Surface aligned ', METHOD)); end
    subplot(4, 6, 6*row);
    pcshow(ptsMatchModel_aligned_disp, 'MarkerSize', 50);
    if row == 1; title(strcat('Model aligned ', METHOD)); end

end

%% helper function that returns the descriptor for a set of local aligned points
% this is copied from getSpacialHistogramDescriptors
function desc = getDesc(pts_local, R)
    
    % descriptor hyperparameters
    NUM_R = 10; % 10
    NUM_THETA = 7; % 7
    NUM_PHI = 14; % 14

    % transform points into spherical coordinates
    r_spheric = vecnorm(pts_local, 2, 2);
    theta_spheric = acos(pts_local(:, 3) ./ r_spheric);
    phi_spheric = atan2(pts_local(:, 2), pts_local(:, 2));
    pts_spheric = [r_spheric, theta_spheric, phi_spheric];

    r_equi = 0:R^3/NUM_R:R^3;
    r_bins = nthroot(r_equi, 3);
    phi_bins = -pi:2*pi/NUM_PHI:pi;
    theta_bins = 0:pi/NUM_THETA:pi;

    % get trivariate histogram
    [counts, ~, ~, ~] = histcn(pts_spheric, r_bins, theta_bins, phi_bins);
    
    % flatten 3D histogram to 1D
    desc = reshape(counts, [], 1);
end

%% helper function that returns the distance between two desciptors
% to the matching that is used in the main script

function dist = getDescDist(desc1, desc2)
    % specify parameters
    change_metric = 0.45;
    unnorm = 2;
    avg_desc_len = 1600; % this is a good estimate obtained by running the main script
    
    % apply unnoramlization
    desc1 = [desc1; unnorm*avg_desc_len];
    desc2 = [desc2; unnorm*avg_desc_len];
    
    % raise to power
    desc1 = desc1.^change_metric;
    desc2 = desc2.^change_metric;    
        
    % normalize
    desc1 = desc1 / norm(desc1, 2);
    desc2 = desc2 / norm(desc2, 2);
    
    % get L1-distance
    dist = vecnorm(desc1 - desc2, 1);    
end

%% helper function: check whether alignment succeeded, given two 
% rotation matrices
% returns true or false based on difference between rotation matrices
% threshhold for success should be between 0.05 and 0.5, depending on strictness. 
% random alignment results in 1.5 < dist < 3 
function dist = checkAlignment(R1, R2)
    dist = norm(R1*R2'-eye(3), 'fro');
    
    %dist = norm(R1*R2'-eye(3), 2);
end

%% helper function: random-uniformly sample point cloud
% custom function, samples the intersection of spherical region for model
% and the cuboid region for the surface. Samples for pcRect are expanded by
% the margin, an samples for pCsphere are shrunk. 
function sample_pts = pcIntersectionSamples(pcRect, pcSphere, d, margin)
    % calculate num_pts to sample based on size of pointcloud and d
    rangeX = pcRect.XLimits(2) - pcRect.XLimits(1) + 2*margin;
    rangeY = pcRect.YLimits(2) - pcRect.YLimits(1) + 2*margin;
    rangeZ = pcRect.ZLimits(2) - pcRect.ZLimits(1) + 2*margin;
    num_pts = round((rangeX * rangeY * rangeZ) / (d^3));
          
    % sample enough random uniformly distributed numbers
    sample_pts = rand(num_pts, 3);
    
    % scale numbers so that they fit into the correct range
    sample_pts = sample_pts .* [rangeX, rangeY, rangeZ];
    sample_pts = sample_pts + ...
        [pcRect.XLimits(1), pcRect.YLimits(1), pcRect.ZLimits(1)] - margin;
    
    % now remove points from region that is not sampled in the Model
    % now get center and only use points that are within "reach"
    center = [(pcSphere.XLimits(2) + pcSphere.XLimits(1))/2, ...
              (pcSphere.YLimits(2) + pcSphere.YLimits(1))/2, ...
              (pcSphere.ZLimits(2) + pcSphere.ZLimits(1))/2];
          
    pts_rel = sample_pts - center;
    dists = vecnorm(pts_rel, 2, 2);
    mask = dists < (rangeX/2 - margin); % note that margin > 0
    
    % based on the mask, return the relevant sample points
    mask = cat(2, mask, mask, mask);
    sample_pts = reshape(sample_pts(mask), [], 3);  
end