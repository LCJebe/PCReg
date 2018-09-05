%% run GetSphericalDescriptors first
close all

% options. If choisng RANSAC, then it's INLIERS_ONLY always!
RANSAC = true;
INLIERS_ONLY = true; % false means: no inliers, only wrong matches. 
MAX_MATCHES = 40;
ALIGN = true;

% figure settings
screensize = get( 0, 'Screensize' );
if ~ALIGN
    figpos = [screensize(3)/3, 75, screensize(3)/3, screensize(4)-150];
else 
    figpos = [screensize(3)/6, 75, 2*screensize(3)/3, screensize(4)-150];
end

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

% for better display: generate the 8 points that span the cube defined by R
add_pts = descOptM.R*[-1,-1,-1;...
                     -1,-1, 1;...
                     -1, 1,-1;...
                     -1, 1, 1;...
                      1,-1,-1;...
                      1,-1, 1;...
                      1, 1,-1;...
                      1, 1, 1];
                  
                  
% visualize random matches
matches_vis = randperm(size(matchesModel, 1));
% remove inliers from matches_vis
matches_vis = setdiff(matches_vis,inlier_idx);
matches_vis = matches_vis(1:num_matches);

for i = 1:num_matches
    if ~RANSAC
        if INLIERS_ONLY
            coordSurface = featSurface(matchesModel(inlier_idx(i), 1), :);
            coordModel = featModel(matchesModel(inlier_idx(i), 2), :);
        else           
            coordSurface = featSurface(matchesModel(matches_vis(i), 1), :);
            coordModel = featModel(matchesModel(matches_vis(i), 2), :);
        end 
    else
        coordSurface = featSurface(matchesModel(inlierPtIdx(i), 1), :);
        coordModel = featModel(matchesModel(inlierPtIdx(i), 2), :);
    end
    
    % get distance between coords (used to determine inliers)
    % transform Surface first, if using RANSAC
    TF = true;
    if RANSAC && TF
        coordModel_aligned = [coordModel, 1]*T;
        coordModel_aligned = coordModel_aligned(1:3);
        coord_dist = norm(coordSurface - coordModel_aligned);
    else
        coord_dist = norm(coordSurface - coordModel);
    end

    
    % get local points for each support region
    ptsMatchSurface = getLocalPoints(ptsSurface, descOptS.R, coordSurface, descOptS.min_pts, descOptS.max_pts);
    ptsMatchModel = getLocalPoints(ptsModel, descOptM.R, coordModel, descOptM.min_pts, descOptM.max_pts);
    
    % optional: align to local reference frame
    if ALIGN
        [ptsMatchSurface_aligned, ~] = AlignPoints_KNN_v2(ptsMatchSurface);
        [ptsMatchModel_aligned, ~] = AlignPoints_KNN_v2(ptsMatchModel);
    end
    
    % add frame points for display
    ptsMatchSurface = [ptsMatchSurface; add_pts];
    ptsMatchModel = [ptsMatchModel; add_pts];
    if ALIGN
        ptsMatchSurface_aligned = [ptsMatchSurface_aligned; add_pts];
        ptsMatchModel_aligned = [ptsMatchModel_aligned; add_pts];
    end
    
    
    % visualize
    if ~ ALIGN
        row = mod(i-1, 4)+1;
        if row ==1
           fig_h = figure();
           set(fig_h,'Position',figpos)
        end
        subplot(4, 2, 2*row-1);
        pcshow(ptsMatchSurface, 'MarkerSize', 50);
        xlabel(sprintf('CoordDist: %0.1f', coord_dist));
        if row == 1; title('Surface'); end
        subplot(4, 2, 2*row);
        pcshow(ptsMatchModel, 'MarkerSize', 50);
        if row == 1; title('Model'); end
    else
        row = mod(i-1, 4)+1;
        if row ==1
           fig_h = figure();
           set(fig_h,'Position',figpos)
        end
        subplot(4, 4, 4*row-3);
        pcshow(ptsMatchSurface, 'MarkerSize', 50);
        xlabel(sprintf('CoordDist: %0.1f', coord_dist));
        if row == 1; title('Surface'); end
        subplot(4, 4, 4*row-2);
        pcshow(ptsMatchModel, 'MarkerSize', 50);
        if row == 1; title('Model'); end
        subplot(4, 4, 4*row-1);
        pcshow(ptsMatchSurface_aligned, 'MarkerSize', 50);
        xlabel(sprintf('%d pts', size(ptsMatchSurface_aligned, 1)));
        if row == 1; title('Surface aligned'); end
        subplot(4, 4, 4*row);
        pcshow(ptsMatchModel_aligned, 'MarkerSize', 50);
        xlabel(sprintf('%d pts', size(ptsMatchModel_aligned, 1)));
        if row == 1; title('Model aligned'); end
    end
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