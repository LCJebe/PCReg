%% run getInliersRANSAC first to load necessary variables into workspace
% or operate on the aligned case!

path = 'Data/PointClouds/';

% options 
RANSAC = true;

%% create one pointcloud that contains matched point pairs post-alignment
if ~ RANSAC
    inliers_vis = inlier_idx;
else
    inliers_vis = inlierPtIdx;
end

% set colors
color1 = repmat([255, 0, 0], size(inliers_vis, 1), 1); % RED: Model
color2 = repmat([0, 255, 0], size(inliers_vis, 1), 1); % GREEN: Surface
color = [color1; color2];

if ~ RANSAC
    pts = single([featModel(matchesModel(inlier_idx, 2), :); ...
                    featSurface(matchesModel(inlier_idx, 1), :)]);
    pcInliers = pointCloud(pts, 'Color', color);
else
    pts = single([pts1_aligned(inliers_vis, :); pts2(inliers_vis, :)]);
    pcInliers = pointCloud(pts, 'Color', color);
end


%% save Inliers
save_file = strcat(path, 'inliers_aligned.pcd');
pcwrite(pcInliers, save_file);

%% DEBUG: Save all matches
% pts = single([pts1_aligned; pts2]);
% 
% color1 = repmat([255, 0, 0], size(pts, 1)/2, 1); % RED: Model
% color2 = repmat([0, 255, 0], size(pts, 1)/2, 1); % GREEN: Surface
% color = [color1; color2];
% 
% pcInliers = pointCloud(pts, 'Color', color);
% 
% save_file = strcat(path, 'matches_aligned.pcd');
% pcwrite(pcInliers, save_file);

%% DEBUG: Save matches after manual alignment with align2

% pts = single([pts1; pts2_postalign]);
% 
% color1 = repmat([255, 0, 0], size(pts, 1)/2, 1); % RED: Model
% color2 = repmat([0, 255, 0], size(pts, 1)/2, 1); % GREEN: Surface
% color = [color1; color2];
% 
% pcInliers = pointCloud(pts, 'Color', color);
% 
% save_file = strcat(path, 'matches_postalign_DELETE.pcd');
% pcwrite(pcInliers, save_file);