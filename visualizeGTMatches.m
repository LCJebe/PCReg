close all

%% load model and aligned surface
path = 'Data/PointClouds/';

pcSurface = pcread(strcat(path, 'Surface_DS3_alignedM.pcd'));
pcModel = pcread(strcat(path, 'GoodCropSpherical_smaller.pcd'));

%% options
MAX_MATCHES = 16;

%% sample points (only intersection of surface and model sampling area)
d = 1;
margin = 3.5;
sample_pts = pcIntersectionSamples(pcSurface, pcModel, d, margin);

%% select points, where surface and model are not empty
R = 3.5;
min_pts = 500;
max_pts = 8000;
thVar = [3, 1.5];

%% determine the valid points and sample only those

valid_mask = false(size(sample_pts, 1), 1);
ptsSurface = pcSurface.Location;
ptsModel = pcModel.Location;

tic
parfor i = 1:size(sample_pts, 1)
    c = sample_pts(i, :);
    [~, distsS] = getLocalPoints(ptsSurface, R, c, min_pts, max_pts);
    [~, distsM] = getLocalPoints(ptsModel, R, c, min_pts, max_pts);
    if ~isempty(distsS) && ~isempty(distsM)
        valid_mask(i) = 1;
    end
end
% retain only valid sample_pts
valid_mask = cat(2, valid_mask, valid_mask, valid_mask);
sample_pts = reshape(sample_pts(valid_mask), [], 3);

valid_pts = [];
valid_numpts = [];
valid_vars = [];
desc_dists = [];
desc_dists2 = [];
desc_dists3 = [];
alignment_meas = [];
alignment_meas2 = [];
alignment_meas3 = [];
centroid_dists2 = [];


for i = 1:size(sample_pts ,1)
    c = sample_pts(i, :);
    % for Surface
    [local_ptsS, ~] = getLocalPoints(pcSurface.Location, R, c, min_pts, max_pts);
    % for Model
    [local_ptsM, ~] = getLocalPoints(pcModel.Location, R, c, min_pts, max_pts);

    % check for variance ratio constraint
    [~, ~, varsS] = pca(local_ptsS, 'Algorithm', 'eig');        
    [~, ~, varsM] = pca(local_ptsM, 'Algorithm', 'eig');
    ratiosS = [varsS(1) / varsS(2), varsS(2) / varsS(3)];
    ratiosM = [varsM(1) / varsM(2), varsM(2) / varsM(3)];
    
    % check if variance constraints are met
    if all((ratiosS - thVar) > 0) && all((ratiosM - thVar) > 0)
        valid_pts = [valid_pts; c];
        valid_numpts = [valid_numpts; size(local_ptsS, 1), size(local_ptsM, 1)];
        valid_vars = [valid_vars; ratiosS, ratiosM];

        % align using both methods and save descriptor distance
        [local_ptsS_aligned, tfS] = AlignPoints(local_ptsS);
        [local_ptsM_aligned, tfM] = AlignPoints(local_ptsM);

        [local_ptsS_aligned2, tfS2, cS2] = AlignPoints_v2(local_ptsS);
        [local_ptsM_aligned2, tfM2, cM2] = AlignPoints_v2(local_ptsM);
        
        [local_ptsS_aligned3, tfS3, cS3] = AlignPoints_KNN(local_ptsS);
        [local_ptsM_aligned3, tfM3, cM3] = AlignPoints_KNN(local_ptsM);

        if ~isempty(local_ptsS_aligned2) && ~isempty(local_ptsM_aligned2)

            % get descriptors for both alignments
            dS = getDesc(local_ptsS_aligned, R);
            dM = getDesc(local_ptsM_aligned, R);

            dS2 = getDesc(local_ptsS_aligned2, R);
            dM2 = getDesc(local_ptsM_aligned2, R);
            
            dS3 = getDesc(local_ptsS_aligned3, R);
            dM3 = getDesc(local_ptsM_aligned3, R);

            % get descriptor distance for both alignments
            dist = getDescDist(dS, dM);
            dist2 = getDescDist(dS2, dM2);
            dist3 = getDescDist(dS3, dM3);

            % save in complete variables
            desc_dists = [desc_dists, dist];
            desc_dists2 = [desc_dists2, dist2];
            desc_dists3 = [desc_dists3, dist3];

            % get measurement for how good the alignment is
            alignment_meas = [alignment_meas, checkAlignment(tfS, tfM)];
            alignment_meas2 = [alignment_meas2, checkAlignment(tfS2, tfM2)];
            alignment_meas3 = [alignment_meas3, checkAlignment(tfS3, tfM3)];
            
            % save centroid distance
            centroid_dists2 = [centroid_dists2, norm(cS2-cM2)];
        else
            % alignment_v2 failed, but we still need to fill the
            % vectors with dummies so that indexing doesn't break
            % "nan" can be ignored with nanmedian and nanmean later!!
            desc_dists = [desc_dists, nan];
            desc_dists2 = [desc_dists2, nan];
            desc_dists3 = [desc_dists3, nan];
            alignment_meas = [alignment_meas, nan];
            alignment_meas2 = [alignment_meas2, nan];
            alignment_meas3 = [alignment_meas3, nan];
            centroid_dists2 = [centroid_dists2, nan];
        end
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
mean_dist = nanmean(desc_dists);
mean_dist2 = nanmean(desc_dists2);
mean_dist3 = nanmean(desc_dists3);
fprintf('Median (Mean) Distance with old alignment: %0.2f (%0.2f)\n', median_dist, mean_dist);
fprintf('Median (Mean) Distance with local alignment: %0.2f (%0.2f)\n', median_dist2, mean_dist2);
fprintf('Median (Mean) Distance with knn alignment: %0.2f (%0.2f)\n', median_dist3, mean_dist3);

%% and more cool statistics: num successes 
successThresh = 0.5;
num_success = sum(alignment_meas < successThresh);
num_success2 = sum(alignment_meas2 < successThresh);
num_success3 = sum(alignment_meas3 < successThresh);
perc_success = 100*num_success/num_valid;
perc_success2 = 100*num_success2/num_valid;
perc_success3 = 100*num_success3/num_valid;
fprintf('Old: %d (%0.1f %%) successes, Local: %d (%0.1f %%) successes, KNN: %d (%0.1f %%) successes\n', ...
            num_success, perc_success, num_success2, perc_success2, num_success3, perc_success3);
        
PLOT = true;
if PLOT
    thresh = 0:0.01:1;
    perc_succ = zeros(length(thresh), 1);
    perc_succ2 = zeros(length(thresh), 1);
    perc_succ3 = zeros(length(thresh), 1);
    for i = 1:length(thresh)
        th = thresh(i);
        num_success = sum(alignment_meas < th);
        num_success2 = sum(alignment_meas2 < th);
        num_success3 = sum(alignment_meas3 < th);
        perc_succ(i) = 100*num_success/num_valid;
        perc_succ2(i) = 100*num_success2/num_valid;
        perc_succ3(i) = 100*num_success3/num_valid;
    end
    plot(thresh, perc_succ);
    hold all
    plot(thresh, perc_succ2);
    plot(thresh, perc_succ3);
    legend('Global PCA', 'Local PCA', 'KNN PCA');
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

    [ptsMatchSurface_aligned2, ~, cS] = AlignPoints_v2(ptsMatchSurface);
    [ptsMatchModel_aligned2, ~, cM] = AlignPoints_v2(ptsMatchModel);

    
    % add frame points for display
    ptsMatchSurface = [ptsMatchSurface; add_pts];
    ptsMatchModel = [ptsMatchModel; add_pts];
    ptsMatchSurface_aligned = [ptsMatchSurface_aligned; add_pts];
    ptsMatchModel_aligned = [ptsMatchModel_aligned; add_pts];
    ptsMatchSurface_aligned2 = [ptsMatchSurface_aligned2; add_pts];
    ptsMatchModel_aligned2 = [ptsMatchModel_aligned2; add_pts];
    
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
    pcshow(ptsMatchSurface_aligned2, 'MarkerSize', 50);
    xlabel(sprintf('align: %0.2f', alignment_meas2(idx)));
    if row == 1; title('Surface aligned v2'); end
    subplot(4, 6, 6*row);
    pcshow(ptsMatchModel_aligned2, 'MarkerSize', 50);
    if row == 1; title('Model aligned v2'); end

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