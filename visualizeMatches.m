%% run GetSphericalDescriptors first
close all

% options. If choisng RANSAC, then it's INLIERS_ONLY always!
RANSAC = false;
INLIERS_ONLY = false;
MAX_MATCHES = 40;

% figure settings
screensize = get( 0, 'Screensize' );
figpos = [screensize(3)/3, 100, screensize(3)/3, screensize(4)-200];

% get coordinates from Model and Surface of Match
if ~RANSAC
    if INLIERS_ONLY
        num_matches = inliers1;
    else
        num_matches = size(matchesModel, 1);
    end
else
    num_matches = size(inlierPtIdx, 1);
end

num_matches = min([num_matches, MAX_MATCHES]);

for i = 1:num_matches
    if ~RANSAC
        if INLIERS_ONLY
            coordSurface = featSurface(matchesModel(inlier_idx(i), 1), :);
            coordModel = featModel(matchesModel(inlier_idx(i), 2), :);
        else
            matches_vis = randperm(size(matchesModel, 1),num_matches);
            coordSurface = featSurface(matchesModel(matches_vis(i), 1), :);
            coordModel = featModel(matchesModel(matches_vis(i), 2), :);
        end
    else
        coordSurface = featSurface(matchesModel(inlierPtIdx(i), 1), :);
        coordModel = featModel(matchesModel(inlierPtIdx(i), 2), :);
    end

    % get local points for each support region
    ptsMatchSurface = getLocalPoints(ptsSurface, R, coordSurface, min_pts);
    ptsMatchModel = getLocalPoints(ptsModel, R, coordModel, min_pts);

    % visualize
    row = mod(i-1, 4)+1;
    if row ==1
       fig_h = figure();
       set(fig_h,'Position',figpos)
    end
    subplot(4, 2, 2*row-1);
    pcshow(pointCloud(ptsMatchSurface), 'MarkerSize', 50);
    if row == 1; title('Surface'); end
    subplot(4, 2, 2*row);
    pcshow(pointCloud(ptsMatchModel), 'MarkerSize', 50);
    if row == 1; title('Model'); end
end

%% DEBUG: get descriptors for those and see if they are actually close...

% [~, desc1, ~] = getMomentDescriptors(ptsMatchSurface, [0,0,0], min_pts, R, thVar, ALIGN_POINTS, CENTER, k);
% [~, desc2, ~] = getMomentDescriptors(ptsMatchModel, [0,0,0], min_pts, R, thVar, ALIGN_POINTS, CENTER, k);
% 
% % normalize
% desc1N = desc1/norm(desc1);
% desc2N = desc2/norm(desc2);
% 
% vecnorm(desc1N-desc2N, 2, 2)
% 
% %% check that distance for other matches (even if they're not inliers)
% 
% descMatchSurface = descSurface(matchesModel(:, 1), :);
% descMatchModel = descModel(matchesModel(:, 2), :);
% 
% descNormsSurface = vecnorm(descMatchSurface, 2, 2);
% descNormsModel = vecnorm(descMatchModel, 2, 2);
% 
% descMatchSurface = descMatchSurface./descNormsSurface;
% descMatchModel = descMatchModel./descNormsModel;
% 
% descMatchDistances = vecnorm(descMatchSurface-descMatchModel, 2, 2);
% 
% %% compare with random distances of all descriptors
% 
% descMatchSurface = descSurface(randperm(100), :);
% descMatchModel = descModel(randperm(100), :);
% 
% descNormsSurface = vecnorm(descMatchSurface, 2, 2);
% descNormsModel = vecnorm(descMatchModel, 2, 2);
% 
% descMatchSurface = descMatchSurface./descNormsSurface;
% descMatchModel = descMatchModel./descNormsModel;
% 
% descMatchDistances = vecnorm(descMatchSurface-descMatchModel, 2, 2);